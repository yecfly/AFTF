// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>



// TODO:  �ڴ˴����ó�����Ҫ������ͷ�ļ�
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