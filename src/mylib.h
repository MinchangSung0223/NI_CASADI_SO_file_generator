#ifndef _MYLIB_H_
#define _MYLIB_H_

#ifdef _WIN32
  #ifdef MYLIB_EXPORTS
    #define MYLIB_API __declspec(dllexport)
  #else
    #define MYLIB_API __declspec(dllimport)
  #endif
#else
  #define MYLIB_API
#endif

#include <Eigen/Dense>  // Eigen 헤더 포함

// 함수 선언: VectorXd의 합을 계산하는 함수
MYLIB_API double vector_sum(const Eigen::VectorXd& vec);

#endif
