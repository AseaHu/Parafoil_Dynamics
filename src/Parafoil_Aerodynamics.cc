#include "Parafoil_Aerodynamics.h"


Parafoil_Aerodynamics::Parafoil_Aerodynamics(){}


/**
 * @brief   翼伞系统初始化：读入配置参数，并计算中间变量
 */
void Parafoil_Aerodynamics::SystemInitialization(){
	ImportPhysicalParameters();
    ImportAerodynamicalParameters();
	ImportInitialState();
	CalConstantParameters();
}


/**
 * @brief   利用读入的数据，计算常值参数
 */
void Parafoil_Aerodynamics::CalConstantParameters(){

	//必要参数
	M = m_p + m_b;
	AR = b / c;

	double tmp0 = M_PI / 4.0;
	k_A = 0.848 * tmp0;
	k_B = 0.339 * tmp0;
	k_C = tmp0 * AR / (1.0 + AR);

	tmp0 = AR / (1.0 + AR);
	k_A_ = 0.055 * tmp0;
	k_B_ = 0.0308 * tmp0;
	k_C_ = 0.0555;
	
	//X_cp = Vector3d( 0.0, 0.0, -5.59 );			//系统质心到翼伞质心矢量在翼伞体坐标系下的表示。在原论文中，该量不明确
	//X_cb = Vector3d( 0.0, 0.0, 0.5 );				//系统质心到载荷质心矢量在载荷体坐标系下的表示。在原论文中，该量不明确
	X_cp = Vector3d( 0.0, 0.0, -5.4 );				//系统质心到翼伞质心矢量在翼伞体坐标系下的表示。在原论文中，该量不明确
	X_cb = Vector3d( 0.0, 0.0, 0.6 );				//系统质心到载荷质心矢量在载荷体坐标系下的表示。在原论文中，该量不明确
	R_cp = SkewSymmetric(X_cp);
	R_cb = SkewSymmetric(X_cb);


	//计算系统惯量矩：翼伞惯量矩I_p，载荷惯量矩I_b
	Matrix3d tmp1 = Matrix3d::Zero();

	///*
	tmp1(0, 0) = b*b + t*t;
	tmp1(1, 1) = c*c + t*t;
	tmp1(2, 2) = b*b + c*c;
	I_p = m_p / 12.0 * tmp1;

	tmp1(0, 0) = y_b*y_b + z_b*z_b;
	tmp1(1, 1) = x_b*x_b + z_b*z_b;
	tmp1(2, 2) = y_b*y_b + x_b*x_b;
	I_b = m_b / 12.0 * tmp1;
	//*/
	/*
	tmp1(0, 0) = m_p * X_cp.z() * X_cp.z();
	tmp1(1, 1) = m_p * X_cp.z() * X_cp.z();
	tmp1(2, 2) = m_p / 12.0 * (b*b + c*c);
	I_p = tmp1;

	tmp1(0, 0) = m_b * X_cb.z() * X_cb.z();
	tmp1(1, 1) = m_b * X_cb.z() * X_cb.z();
	tmp1(2, 2) = m_b / 12.0 * (y_b*y_b + x_b*x_b);
	I_b = tmp1;
	*/

	I = I_p + I_b;


	//计算初始状态
	phi = phi / 180.0 * M_PI;						//读入的角度单位均为：degree
	theta = theta / 180.0 * M_PI;					//读入的角度单位均为：degree
	psi = psi / 180.0 * M_PI;						//读入的角度单位均为：degree
	gliding_angle = gliding_angle / 180.0 * M_PI;	//读入的角度单位均为：degree
	omega = omega / 180.0 * M_PI;					//读入的角度单位均为：degree
	Matrix3d R_v2b = AngleAxisd(gliding_angle - theta, Vector3d(0, 1, 0)).toRotationMatrix();
	V_c = R_v2b * Vector3d(velocity, 0.0, 0.0);		//初始速度
}


/**
 * @brief   根据当前海拔高度，计算空气密度
 * @param	altitude
 */
void Parafoil_Aerodynamics::CalAirDensity(double altitude){
	
	if(altitude < 11000){
		double tmp = 1 - altitude / 44330.0;
		tho = tho0 * pow(tmp, 4.256);
	}
	else if(altitude < 20000){
		double tmp = - ( altitude - 11000.0) / R_e;
		tho = 0.3639 * exp(tmp);
	}
}


/**
 * @brief   系统平移、姿态运动方程，输出位置、姿态角的导数
 * @param  	V
 * @param  	F
 * @param  	omega
 */
Matrix<double, 6, 1> Parafoil_Aerodynamics::MotionFunction(Matrix<double, 6, 1> &state){
	Matrix<double, 6, 1> dot_f;

	double phi = state(3);
	double theta = state(4);
	double psi = state(5);

	Vector3d angle = state.block(3,0,3,1);
	Matrix3d T = EulerAngle2RotationMatrix(angle);		//Ti表示将参考系下的向量转换为本体系下的表示

	Matrix3d F;
	F<<	1,	sin(phi) * tan(theta), 	cos(phi) * tan(theta), 
		0,	cos(phi),				-sin(phi),
		0,	sin(phi) / cos(theta),	cos(phi) / cos(theta);

	dot_f.block(0,0,3,1) = T.transpose() * V_c;
	dot_f.block(3,0,3,1) = F * omega;

	return dot_f;
}

/**
 * @brief   系统加速度、角加速度方程，输出速度、角速度的导数
 * @param  	state
 */
Matrix<double, 6, 1> Parafoil_Aerodynamics::AccelerationFunction(Matrix<double, 6, 1> &state){
	Vector3d V_c = state.block(0,0,3,1);
	Vector3d omega = state.block(3,0,3,1);

	Matrix3d Omega = SkewSymmetric(omega);
	V_p = V_c + Omega * X_cp;
	V_b = V_c + Omega * X_cb;
	

	//更新作用于翼伞和载荷的气动力
	Vector3d ang(phi, theta, psi);
	Matrix3d T = EulerAngle2RotationMatrix(ang);	//表示将参考系下的向量转换为本体系下的表示
	Vector3d vp = V_p - T * V_wind;
	Vector3d vb = V_b - T * V_wind;
	Vector3d vc = V_c - T * V_wind;
	Matrix3d Xi_p = SkewSymmetric(vp);

	alpha_p = atan2(vp.z(), vp.x());		//翼伞的攻角
	C_L = C_L_0 + C_L_alpha * alpha_p;
	C_D_p = C_D_0 + C_D_alpha * alpha_p * alpha_p;
	F_A_p = tho * S_p * vp.norm() * ( C_L * Vector3d(vp.z(), 0.0, -vp.x()) - C_D_p * vp) / 2.0;

	
	alpha_b = atan2(vb.z(), vb.x());		//载荷的攻角
	C_D_b = C_D_0 + C_D_alpha * alpha_b * alpha_b;
	F_A_b = - tho * S_b * vb.norm() * C_D_b * vb / 2.0;


	//更新造成系统滚转、俯仰、偏航运动的气动力矩 M_A_p
	double Vp = vp.norm();
	double tmpx = tho * S_p * Vp * Vp * C_L / 2.0;
	double x = C_l_p * b * b * omega.x() / Vp / 2.0 + C_l_phi * b * phi;
	double y = C_m_q * c * c * omega.y() / Vp / 2.0 + C_m_0 * c + C_m_alpha * c * alpha_p;
	double z = C_n_r * b * b * omega.z() / Vp / 2.0;

	M_A_p = tmpx * Vector3d( x, y, z );


	//更新由襟翼偏转造成的气动力增益矩阵S_F_ap
	double delta_s = min(delta_left, delta_right);
	double delta_a = delta_left - delta_right;

	double tmpd = 1;
	if( delta_a < 0)
		tmpd = -1;
	Matrix<double, 3, 2> Sc;
	Sc<<	tmpd * ( C_L_da * vp.z() - C_D_da * vp.x() ),			C_L_ds * vp.z() - C_D_ds * vp.x(),
			-tmpd * C_D_da * vp.y(),								-C_D_ds * vp.y(),
			tmpd * ( C_L_da * vp.x() - C_D_da * vp.z() ),			C_L_ds * vp.x() - C_D_ds * vp.z();
	S_F_ap = tho * S_p * Vp / 2.0 * Sc;


	//更新由襟翼偏转造成的气动力矩增益矩阵S_M_ap
	Matrix<double, 3, 2> Mc;
	Mc<<	C_l_da * b / t,		0,
			0,					0,
			C_n_da * b / t,		0;
	S_M_ap = tho * S_p * Vp * Vp / 2.0 * Mc;


	//更新翼伞系统动力学方程（状态转移矩阵和控制矩阵等）
	Matrix<double, 6, 1> E = Matrix<double, 6, 1>::Zero();
	E.block(0,0,3,1) = W_g + F_A_p + F_A_b - Omega * M_F * vp - (M * Matrix3d::Identity() + M_F) * Omega * vc;
	E.block(3,0,3,1) = M_A_p - Xi_p * M_F * vp + R_cp * F_A_p - R_cp * Omega * M_F * vp + R_cb * F_A_b - Omega * (I + I_F) * omega;
	A = D.inverse() * E;

	S.block(0,0,3,2) = S_F_ap;
	S.block(3,0,3,2) = S_M_ap + R_cp * S_F_ap;
	B = D.inverse() * S;

	Matrix<double, 6, 1> dot_f = A + B * Vector2d( delta_a, delta_s );

	return dot_f;
}


/**
 * @brief   根据翼伞参数及输入的控制参数(襟翼偏转角度)，建立翼伞动力学模型并模拟翼伞飞行轨迹，输出翼伞状态
 * @param 	d_left		左侧襟翼偏转量（°）
 * @param	d_right		右侧襟翼偏转量（°）
 */
PVA Parafoil_Aerodynamics::Sim_With_SingleStep(double d_left, double d_right, const Vector3d & v_wind){
	PVA pva;
	delta_left = d_left / 180.0 * M_PI;
	delta_right = d_right / 180.0 * M_PI;
	V_wind = v_wind;

	UpdateCurrentStates();

	pva.Position = Position;
	pva.Velocity = V_w;

	return pva;
}


/**
 * @brief   根据翼伞参数及导入的控制参数(包括仿真时间、控制量等)，建立翼伞动力学模型并模拟翼伞飞行轨迹
 */
void Parafoil_Aerodynamics::Sim_With_DesignedInputs(){

	ImportControlCommands();

	int k = 0;
	int times = int(SimTime / TimeStep);
	int node = int(Commands[k].x() / TimeStep);

	double angle = gliding_angle * M_PI / 180.0;
	Positions.push_back(Position);
	Velocities.push_back(velocity * Vector3d(cos(angle), 0.0, -sin(angle)));
	EulerAngles.push_back(Vector3d(phi, theta, psi));
	AngularRates.push_back(omega);

	cout<<"-----------------------Control Command Sequences-----------------------"<<endl;
	cout<<setw(10)<<"Time"<<setw(30)<<"Left flap delection/deg"<<setw(30)<<"Right flap delection/deg"<<endl;
	for(int i = 0; i < times; i ++){
		if(i == node){
			delta_left = Commands[k].y() / 180.0 * M_PI;
			delta_right = Commands[k].z() / 180.0 * M_PI;
			cout<<setw(10)<<Commands[k].x()<<setw(30)<<delta_left*180.0/M_PI<<setw(30)<<delta_right*180.0/M_PI<<endl;

			k++;
			if(k < Commands.size())
				node = int(Commands[k].x() / TimeStep);
		}
		UpdateCurrentStates();
	}
}


/**
 * @brief   根据当前状态，更新中间变量，得到动力学方程中的系统状态转移矩阵和控制矩阵
 */
void Parafoil_Aerodynamics::UpdateCurrentStates(){

	//更新空气密度
	double altitude = -Position.z();
	CalAirDensity(altitude);

	//更新附加质量矩阵 M_F
	Matrix3d tmp = Matrix3d::Zero();
	tmp(0, 0) = k_A * tho * t * t * b * ( 1.0 + 8.0 / 3.0 * a * a * a );
	tmp(1, 1) = k_B * tho * t * t * ( t * t + 2.0 * a * a * ( 1 - t * t ) ) * c;
	tmp(2, 2) = k_C * tho * c * c * b * sqrt( 1.0 + 2.0 * a * a * ( 1.0 - t * t ) );
	M_F = tmp;


	//更新附加质量惯量矩 I_F
	tmp(0, 0) = k_A_ * tho * c * c * b * b * b;
	tmp(1, 1) = k_B_ * tho * c * c * c * c * b * ( 1.0 + M_PI / 6.0 * (1.0 + AR) * AR * a * a * t * t);
	tmp(2, 2) = k_C_ * tho * t * t * b * b * b * ( 1.0 + 8.0 * a * a );
	I_F = tmp;

	//更新系统重力在系统体坐标系下的表示
	W_g = M * g * Vector3d( -sin(theta), sin(phi) * cos(theta), cos(phi) * cos(theta) );


	//更新翼伞系统动力学方程（状态转移矩阵和控制矩阵等）
	D = Matrix<double, 6, 6>::Zero();
	D.block(0,0,3,3) = M * Matrix3d::Identity() + M_F;
    D.block(3,0,3,3) = R_cp * M_F;
    D.block(3,3,3,3) = I + I_F;


	//四阶龙格库塔法求 V_c 和 omega 积分
	Matrix<double, 6, 1> y0;
	y0.block(0,0,3,1) = V_c;
	y0.block(3,0,3,1) = omega;
	RungeKutta_VO(y0, TimeStep);


	//四阶龙格库塔法求 Position 和 EulerAngle 积分
	Matrix<double, 6, 1> y1;
	y1.block(0,0,3,1) = Position;
	y1.block(3,0,3,1) = Vector3d(phi, theta, psi);
	RungeKutta_PE(y1, TimeStep);


	//更新状态：位置、姿态、速度、角速度
	V_c = y0.block(0,0,3,1);
	omega = y0.block(3,0,3,1);
	Position = y1.block(0,0,3,1);
	phi = y1(3);
	theta = y1(4);
	psi = y1(5);

	Vector3d angle = y1.block(3,0,3,1);
	Matrix3d T = EulerAngle2RotationMatrix(angle);		//Ti表示将参考系下的向量转换为本体系下的表示
	V_w = T.transpose() * V_c;							//转换到参考系下的速度

	CurTime += TimeStep;

	//输出
	Positions.push_back(Position);
	Velocities.push_back(V_w);
	EulerAngles.push_back(Vector3d(phi, theta, psi) * 180.0 / M_PI);
	AngularRates.push_back(omega * 180.0 / M_PI);
}


/**
 * @brief   返回三维向量v的反对称矩阵
 * @param	v
 */
Matrix3d SkewSymmetric(Vector3d &v){
	Matrix3d result;
	result<<	0,			-v(2),		v(1),
				v(2),		0,			-v(0),
				-v(1),		v(0),		0;
	return result;
}


/**
 * @brief   返回欧拉角angle所对应的旋转矩阵，3-2-1旋转序列
 * @param	angle
 */
Matrix3d EulerAngle2RotationMatrix(Vector3d &angle){
	Matrix3d R;

	double phi = angle.x();
	double theta = angle.y();
	double psi = angle.z();

	//R表示将参考系下的向量转换为本体系下的表示
	R << 	cos(theta)*cos(psi),									cos(theta)*sin(psi),									-sin(theta),
			sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi),		sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi),		sin(phi)*cos(theta),
			cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi),		cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi),		cos(phi)*cos(theta);
	return R;
}


/**
 * @brief   四阶龙格库塔法，求速度(V_c)和角速度(omega)的积分
 * @param	y
 * @param	step
 */
void Parafoil_Aerodynamics::RungeKutta_VO(Matrix<double, 6, 1> &y, double step){
	Matrix<double, 6, 1> y1, y2, y3, k1, k2, k3, k4;

	k1 = AccelerationFunction(y);

	y1 = y + step / 2.0 * k1;
	k2 = AccelerationFunction(y1);

	y2 = y + step / 2.0 * k2;
	k3 = AccelerationFunction(y2);

	y3 = y + step * k3;
	k4 = AccelerationFunction(y3);

	y += step / 6.0 * ( k1 + 2*k2 + 2*k3 + k4 );
}


/**
 * @brief   四阶龙格库塔法，求位置(Position)和姿态角(Euler Angle)的积分
 * @param	y
 * @param	step
 */
void Parafoil_Aerodynamics::RungeKutta_PE(	Matrix<double, 6, 1> &y, double step){
	Matrix<double, 6, 1> y1, y2, y3, k1, k2, k3, k4;

	k1 = MotionFunction(y);

	y1 = y + step / 2.0 * k1;
	k2 = MotionFunction(y1);

	y2 = y + step / 2.0 * k2;
	k3 = MotionFunction(y2);

	y3 = y + step * k3;
	k4 = MotionFunction(y3);

	y += step / 6.0 * ( k1 + 2*k2 + 2*k3 + k4 );
}


/**
 * @brief   从parameters文件夹中读取aerodncs_parameters.txt文件并赋值给对应的变量
 */
void Parafoil_Aerodynamics::ImportAerodynamicalParameters(){
	cv::FileStorage fs;
    fs.open("../config/param_aerodynamics.yaml", cv::FileStorage::READ);

	fs["C_L_0"]>>C_L_0;
	fs["C_D_0"]>>C_D_0;
	fs["C_l_p"]>>C_l_p;
	fs["C_m_q"]>>C_m_q;
	fs["C_n_r"]>>C_n_r;
	fs["C_L_da"]>>C_L_da;
	fs["C_D_da"]>>C_D_da;
	fs["C_l_da"]>>C_l_da;
	fs["C_L_alpha"]>>C_L_alpha;
	fs["C_D_alpha"]>>C_D_alpha;
	fs["C_l_phi"]>>C_l_phi;
	fs["C_m_0"]>>C_m_0;
	fs["C_m_alpha"]>>C_m_alpha;
	fs["C_L_ds"]>>C_L_ds;
	fs["C_D_ds"]>>C_D_ds;
	fs["C_n_da"]>>C_n_da;

	cout<<"------------------Aerodynamical Parameters of Parafoil-----------------"<<endl;
	cout<<"C_L_0			"<<C_L_0<<endl;
	cout<<"C_D_0			"<<C_D_0<<endl;
	cout<<"C_l_p			"<<C_l_p<<endl;
	cout<<"C_m_q			"<<C_m_q<<endl;
	cout<<"C_n_r			"<<C_n_r<<endl;
	cout<<"C_L_da			"<<C_L_da<<endl;
	cout<<"C_D_da			"<<C_D_da<<endl;
	cout<<"C_l_da			"<<C_l_da<<endl;
	cout<<"C_L_alpha		"<<C_L_alpha<<endl;
	cout<<"C_D_alpha		"<<C_D_alpha<<endl;
	cout<<"C_l_phi			"<<C_l_phi<<endl;
	cout<<"C_m_0			"<<C_m_0<<endl;
	cout<<"C_m_alpha		"<<C_m_alpha<<endl;
	cout<<"C_L_ds			"<<C_L_ds<<endl;
	cout<<"C_D_ds			"<<C_D_ds<<endl;
	cout<<"C_n_da			"<<C_n_da<<endl;
}


/**
 * @brief   从parameters文件夹中读取physical_parameters.txt文件并赋值给对应的变量
 */
void Parafoil_Aerodynamics::ImportPhysicalParameters(){
	cv::FileStorage fs;
    fs.open("../config/param_physical.yaml", cv::FileStorage::READ);

	fs["m_p"]>>m_p;
	fs["m_b"]>>m_b;
	fs["S_p"]>>S_p;
	fs["S_b"]>>S_b;
	fs["a"]>>a;
	fs["b"]>>b;
	fs["t"]>>t;
	fs["c"]>>c;
	fs["x_b"]>>x_b;
	fs["y_b"]>>y_b;
	fs["z_b"]>>z_b;

	cout<<"---------------------Physical Parameters of System---------------------"<<endl;
	cout<<"m_p          		"<<m_p<<endl;
	cout<<"m_b          		"<<m_b<<endl;
	cout<<"S_p          		"<<S_p<<endl;
	cout<<"S_b          		"<<S_b<<endl;
	cout<<"a            		"<<a<<endl;
	cout<<"b            		"<<b<<endl;
	cout<<"t            		"<<t<<endl;
	cout<<"c            		"<<c<<endl;
	cout<<"x_b          		"<<x_b<<endl;
	cout<<"y_b          		"<<y_b<<endl;
	cout<<"z_b          		"<<z_b<<endl;
}


/**
 * @brief   从parameters文件夹中读取Initial_State.txt文件并赋值给对应的变量
 */
void Parafoil_Aerodynamics::ImportInitialState(){
	cv::FileStorage fs;
    fs.open("../config/init_state.yaml", cv::FileStorage::READ);

	cv::FileNode P = fs["Position"];
	P["x"]>>Position(0);
	P["y"]>>Position(1);
	P["z"]>>Position(2);

	fs["Velocity"]>>velocity;
	fs["GlidingAngle"]>>gliding_angle;

	cv::FileNode EulerAngle = fs["EulerAngle"];
	EulerAngle["roll"]>>phi;
	EulerAngle["pitch"]>>theta;
	EulerAngle["yaw"]>>psi;

	cv::FileNode AngularRate = fs["AngularRate"];
	AngularRate["wx"]>>omega(0);
	AngularRate["wy"]>>omega(1);
	AngularRate["wz"]>>omega(2);

	fs["TimeStep"]>>TimeStep;

	cout<<"------------------------Initial State of System------------------------"<<endl;
	cout<<"Position	"<<Position.transpose()<<endl;
	cout<<"velocity		"<<velocity<<endl;
	cout<<"gliding_angle		"<<gliding_angle<<endl;
	cout<<"EulerAngle		"<<Vector3d(phi, theta, psi).transpose()<<endl;
	cout<<"omega			"<<omega.transpose()<<endl;
	cout<<"TimeStep		"<<TimeStep<<endl;


	Positions.clear();
	Velocities.clear();
	EulerAngles.clear();
	AngularRates.clear();
}


/**
 * @brief   从parameters文件夹中读取Control_Commands.txt文件并赋值给对应的变量
 */
void Parafoil_Aerodynamics::ImportControlCommands(){
	string path = "../config/Control_Commands.txt";

    ifstream TXTReader;
	TXTReader.open(path, ios::in);
	if(!TXTReader){
		cerr<<"Control_Commands.txt open error!"<<endl; 
		abort();
	}

	Commands.clear();

	TXTReader>>SimTime;

	Vector3d tmp;
	while (!TXTReader.eof())
	{
		TXTReader>>tmp(0);
		TXTReader>>tmp(1);
		TXTReader>>tmp(2);
		Commands.push_back(tmp);
	}
	Commands.pop_back();
	//cout<<"Commands: "<<Commands.size()<<"		"<<Commands[0].transpose()<<endl;
	TXTReader.close();
}


/**
 * @brief   返回系统当前时间
 */
double Parafoil_Aerodynamics::Get_Time(){
	return CurTime;
}


/**
 * @brief   保存系统的状态（位置、速度、姿态、角速度）
 */
void Parafoil_Aerodynamics::Save_State(){
	ofstream writer("../outputs/states.txt");
	for(int i = 0; i < Positions.size(); i ++){
		double time = double(i) * TimeStep;
		while(EulerAngles[i](2)>180)
			EulerAngles[i](2) -= 360;
		while(EulerAngles[i](2)<-180)
			EulerAngles[i](2) += 360;
		///*
		writer<<setprecision(12)<<time<<",";
		writer<<setprecision(12)<<Positions[i](0)<<","<<Positions[i](1)<<","<<Positions[i](2)<<",";
		writer<<setprecision(12)<<Velocities[i](0)<<","<<Velocities[i](1)<<","<<Velocities[i](2)<<",";
		writer<<setprecision(12)<<EulerAngles[i](0)<<","<<EulerAngles[i](1)<<","<<EulerAngles[i](2)<<",";
		writer<<setprecision(12)<<AngularRates[i](0)<<","<<AngularRates[i](1)<<","<<AngularRates[i](2)<<endl;
		//*/
		/*
		int precision = 9, width = 16;
		writer<<setprecision(precision)<<setw(4)<<time;
		writer<<setprecision(precision)<<setw(width)<<Positions[i](0)<<setw(width)<<Positions[i](1)<<setw(width)<<Positions[i](2);
		writer<<setprecision(precision)<<setw(width)<<Velocities[i](0)<<setw(width)<<Velocities[i](1)<<setw(width)<<Velocities[i](2);
		writer<<setprecision(precision)<<setw(width)<<EulerAngles[i](0)<<setw(width)<<EulerAngles[i](1)<<setw(width)<<EulerAngles[i](2);
		writer<<setprecision(precision)<<setw(width)<<AngularRates[i](0)<<setw(width)<<AngularRates[i](1)<<setw(width)<<AngularRates[i](2)<<endl;
		*/
	}
	writer.close();
}