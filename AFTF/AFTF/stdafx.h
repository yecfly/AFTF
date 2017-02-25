// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO:  在此处引用程序需要的其他头文件
/*session1 essential, needed in EstRXY class*/
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <cv.h>
#include <opencv/cxcore.hpp>
#include <opencv.hpp>
#include <vector>

#include <stdlib.h>
#include <io.h>
#include <string.h>
#include <windows.h>

//#define showProcess
//#define ResizeToHalf

//#define NormalizeInputImages


struct ImageSet {
	std::vector<std::string> dir;//frame directory
	std::vector<cv::Mat> fc;//frame contect
};
#include "AFTFutl.h"
#include "EstRXY.h"