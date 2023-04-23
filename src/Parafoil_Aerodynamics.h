#ifndef _PARAFOIL_AERODYNAMICS_H_
#define _PARAFOIL_AERODYNAMICS_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace Eigen;


/**
 * @brief   位置、速度和姿态
 */
struct PVA{
    Vector3d Position = Vector3d(0, 0, 0);
    Vector3d Velocity = Vector3d(0, 0, 0);
    Vector3d Attitude = Vector3d(0, 0 ,0);
};


/**
 * @brief   建立翼伞类，构建翼伞系统的六自由度动力学模型，提供输入（初始状态量+控制量）/ 输出（系统状态序列）接口
 */
class Parafoil_Aerodynamics{

public:
    Parafoil_Aerodynamics();
    void SystemInitialization();
    void Sim_With_DesignedInputs();
    PVA Sim_With_SingleStep(double d_left, double d_right, const Vector3d & v_wind);
    double Get_Time();
    void Save_State();

private:
    void ImportPhysicalParameters();
    void ImportAerodynamicalParameters();
    void ImportInitialState();
    void CalConstantParameters();
    void CalAirDensity(double altitude);
    Matrix<double, 6, 1> AccelerationFunction(Matrix<double, 6, 1> &state);
    Matrix<double, 6, 1> MotionFunction(Matrix<double, 6, 1> &state);
    void RungeKutta_VO(Matrix<double, 6, 1> &y, double step);
    void RungeKutta_PE(Matrix<double, 6, 1> &y, double step);
    void UpdateCurrentStates();
    void ImportControlCommands();

public:
    /*************************翼伞系统的状态量*************************/
    double phi;                     //系统姿态欧拉角 - 滚转角，rad
    double theta;                   //系统姿态欧拉角 - 俯仰角，rad
    double psi;                     //系统姿态欧拉角 - 偏航角，rad
    double alpha_p;                 //翼伞的攻角，rad
    double alpha_b;                 //载荷的攻角，rad
    double velocity;                //初始速度大小, m/s
    double gliding_angle;           //初始下滑角, 读入的单位：degree
    Vector3d Position;              //系统质心的位置，表示在世界系下, m
    Vector3d V_w;                   //系统质心的速度(表示在世界参考系下), m/s
    Vector3d V_c;                   //系统质心的速度(表示在本体系下), m/s
    Vector3d V_p;                   //翼伞质心的速度(表示在翼伞本体系下), m/s
    Vector3d V_b;                   //载荷质心的速度(表示在载荷本体系下), m/s
    Vector3d omega;                 //翼伞系统在翼伞体坐标系下的角速度，rad/s？注：论文没有给出明确的定义
    

private:
    /***************************翼伞气动参数***************************/
    double C_L_0, C_D_0, C_l_p, C_m_q, C_n_r, C_L_da, C_D_da, C_l_da, C_L_alpha, C_D_alpha, C_l_phi, C_m_0, C_m_alpha, C_L_ds, C_D_ds, C_n_da;
    double C_L, C_D_p, C_D_b;

    /*************翼伞系统物理参数（尺寸、质量、面积、惯量矩等）*************/
    double a;                       //翼伞弧面高度。注：原文未给出确定值，其影响较大
    double b;                       //翼伞展长
    double c;                       //翼伞弦长
    double t;                       //翼伞厚度
    double x_b;                     //载荷长
    double y_b;                     //载荷宽
    double z_b;                     //载荷高
    double m_p;                     //翼伞质量
    double m_b;                     //载荷质量
    double S_p;                     //翼伞面积
    double S_b;                     //载荷横截面积
    double M;                       //系统总质量
    Vector3d X_cp;                  //系统质心到翼伞质心矢量在翼伞体坐标系下的表示。在原论文中，该量不明确
    Vector3d X_cb;                  //系统质心到载荷质心矢量在载荷体坐标系下的表示。在原论文中，该量不明确
    Matrix3d M_F;                   //附加质量矩阵
    Matrix3d I_p;                   //翼伞惯量矩
    Matrix3d I_b;                   //载荷惯量矩
    Matrix3d I_F;                   //附加质量惯量矩
    Matrix3d I;                     //系统总惯量矩 I = I_p + I_b
    Matrix3d R_cp;                  //系统质心到翼伞质心矢量在翼伞体坐标系下的反对称矩阵
    Matrix3d R_cb;                  //系统质心到载荷质心矢量在载荷体坐标系下的反对称矩阵
    Matrix<double, 3, 2> S_F_ap;    //由襟翼偏转造成的气动力增益矩阵
    Matrix<double, 3, 2> S_M_ap;    //由襟翼偏转造成的气动力矩增益矩阵

    double AR;                      //展弦比
    double k_A, k_B, k_C;           //计算附加质量项M_F的中间参数
    double k_A_, k_B_, k_C_;        //计算附加质量惯量矩I_F的中间参数


    /**************************翼伞系统受力***************************/
    Vector3d W_g;                   //系统重力矢量在系统体坐标系下的表示
    Vector3d F_A_p;                 //作用于翼伞的气动力, 表示在翼伞本体系下
    Vector3d F_A_b;                 //作用于载荷的气动力, 表示在载体本体系下
    Vector3d M_A_p;                 //造成系统滚转、俯仰、偏航运动的气动力矩


    /****************************环境参数****************************/
    //常量
    const double tho0 = 1.2250;     //海平面处的空气密度, kg/m^3
    const double g = 9.8;           //重力加速度，m/s^2
    const double R_e = 6371e3;		//地球半径, m

    //变量
    double tho;                     //当前高度下的空气密度, kg/m^3
    Vector3d V_wind;                //风速，表示在世界系下，m/s

    /****************************控制参数****************************/
    double delta_left = 0;          //左侧襟翼偏转量，rad
    double delta_right = 0;         //右侧襟翼偏转量，rad


    /****************************仿真相关****************************/
    vector<Vector3d> Positions;     //位置序列
    vector<Vector3d> Velocities;    //速度序列
    vector<Vector3d> EulerAngles;   //姿态序列
    vector<Vector3d> AngularRates;  //角速度序列
    vector<Vector3d> Commands;      //控制指令：time, delta_left, delta_right
    double SimTime;                 //仿真时长, 单位：s
    double TimeStep;                //时间步长, 单位：s
    double CurTime = 0;             //当前时间, 单位：s

    /****************************其他参数****************************/
    Matrix<double, 6, 1> A;         //动力学方程中的系统转移矩阵
    Matrix<double, 6, 2> B;         //动力学方程中的控制矩阵
    Matrix<double, 6, 6> D;         //系统质量、转动惯量矩阵
    Matrix<double, 6, 2> S;         //襟翼偏转到输出力和力矩的增益矩阵
};


Matrix3d SkewSymmetric(Vector3d &v);
Matrix3d EulerAngle2RotationMatrix(Vector3d &angle);


#endif

