#include <iostream>
#include "Eigen/Dense"    // Eigen头文件，<Eigen/Dense>包含Eigen库里面所有的函数和类
 
int main()
{

//   test1: input a matrix and output it

  Eigen::MatrixXd m(2,2);   // MatrixXd 表示的是动态数组，初始化的时候指定数组的行数和列数
  m(0,0) = 3;               //m(i,j) 表示第i行第j列的值，这里对数组进行初始化
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;     // eigen重载了<<运算符，可以直接输出eigen矩阵的值


//   test2: matrix operation test
 
  Eigen::Matrix2d a;
  a << 1, 2,
       3, 4;
  Eigen::MatrixXd b(2,2); // 定义动态矩阵--初始化动态数组大小
  b << 2, 3,
       1, 4;
  std::cout << "a + b =\n" << a + b << std::endl;  //矩阵乘法
  std::cout << "a - b =\n" << a - b << std::endl;  //矩阵减法
  std::cout << "Doing a += b;" << std::endl;        //矩阵累加
  a += b;
  std::cout << "Now a =\n" << a << std::endl;
  Eigen::Vector3d v(1,2,3);
  Eigen::Vector3d w(1,0,0);
  std::cout << "-v + w - v =\n" << -v + w - v << std::endl;
}
