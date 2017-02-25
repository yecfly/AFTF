// AFTF.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <direct.h>
#include "EstRXY.h"

using namespace cv;
using namespace std;

////
//from tencent coder's
using namespace cv;
float maxError = 0.06f;

void KeypointsToPoints(const std::vector<KeyPoint> &keyPoints, std::vector<Point2f> &points)
{
	for (int i = 0; i < keyPoints.size(); i++)
	{
		points.push_back(keyPoints[i].pt);
	}
}
/*
float CalLength(Point3f p)
{
	return  sqrtf(p.x *p.x + p.y * p.y + p.z *p.z);
}

inline float sign(float t)
{
	if (t >= 0.0f)
		return 1.0f;
	else
		return -1.0f;
}

std::vector<float> CalBarycentric(const std::vector<Point2f> &points, Point2f tf)
{
	Point3f p1 = Point3f(points[0].x, points[0].y, 0.0f);
	Point3f p2 = Point3f(points[1].x, points[1].y, 0.0f);
	Point3f p3 = Point3f(points[2].x, points[2].y, 0.0f);

	Point3f f(tf.x, tf.y, 0.0f);
	Point3f f1 = p1 - f;
	Point3f f2 = p2 - f;
	Point3f f3 = p3 - f;

	Point3f va = (p1 - p2).cross(p1 - p3);
	Point3f va1 = f2.cross(f3);
	Point3f va2 = f3.cross(f1);
	Point3f va3 = f1.cross(f2);

	float a = CalLength(va);
	float a1 = CalLength(va1) / a * sign(va.dot(va1));
	float a2 = CalLength(va2) / a * sign(va.dot(va2));
	float a3 = CalLength(va3) / a * sign(va.dot(va3));

	std::vector<float> w;
	w.push_back(a1);
	w.push_back(a2);
	w.push_back(a3);
	return w;
}

void tCalRotationAndTranslation(const std::vector<KeyPoint> &object0, const std::vector<KeyPoint> &object1, const std::vector< DMatch > &matches,
	float &rotationAngle, float &dis)
{
	std::vector<DMatch> copyMatches(matches);

	std::cout << copyMatches[0].distance << "  " << copyMatches[1].distance << std::endl;
	std::vector<float> thetaValue;
	for (int i = 0; i < 3; i++)
	{
		if (copyMatches[i + 1].distance > maxError)
			break;
		cv::Point2f d0 = (object0[copyMatches[i].queryIdx].pt - object0[copyMatches[i + 1].queryIdx].pt);
		cv::Point2f d1 = (object1[copyMatches[i].trainIdx].pt - object1[copyMatches[i + 1].trainIdx].pt);
		double len = sqrtf(d0.x * d0.x + d0.y * d0.y);
		d0.x /= len;
		d0.y /= len;

		double len1 = sqrtf(d1.x * d1.x + d1.y * d1.y);
		d1.x /= len1;
		d1.y /= len1;
		//std::cout << len << "  " << len1 << std::endl;

		double cosTheta = d0.dot(d1);
		double sinTheta = d0.cross(d1);

		float tTheta = acos(cosTheta) / 3.1415926 * 180.0;
		if (sinTheta < 0.0f)
			tTheta = -tTheta;
		thetaValue.push_back(tTheta);

	}

	float sumTheta = 0.0f;
	for (int i = 0; i < thetaValue.size(); i++)
	{
		sumTheta += thetaValue[i];
	}

	rotationAngle = sumTheta / thetaValue.size();

	cv::Point2f centerPos(320, 180);
	cv::Point2f translation;

	int count = 0;
	for (int i = 0; i < 3; i++)
	{
		if (copyMatches[i].distance > maxError)
			continue;

		cv::Point2f dd0 = object0[copyMatches[i].queryIdx].pt - centerPos;

		float angleInRadian = rotationAngle / 180.0f * 3.1415926f;
		cv::Point2f nDD0;
		nDD0.x = dd0.x * cosf(angleInRadian) - dd0.y * sinf(angleInRadian);
		nDD0.y = dd0.x * sinf(angleInRadian) + dd0.y * cosf(angleInRadian);

		cv::Point2f newDD0 = centerPos + nDD0;
		cv::Point2f tv = newDD0 - object1[copyMatches[i].trainIdx].pt;
		//std::cout << "Sub_translation:  " << tv << std::endl;
		translation += tv;
		count++;
	}
	translation.x /= count;
	translation.y /= count;

	std::cout << "translation:  " << translation << std::endl;

	dis = 0.0f;
}

void CalMatchesThroughOpticalFlk(const Mat &img1, const Mat &img2,
	const std::vector<KeyPoint> &left_keypoints, const std::vector<KeyPoint> &right_keypoints,
	std::vector<DMatch> &matches)
{
	std::vector<Point2f> left_points;
	KeypointsToPoints(left_keypoints, left_points);

	std::vector<Point2f> right_points(left_points.size());

	std::vector<uchar> vstatus;
	std::vector<float> verror;

	cv::calcOpticalFlowPyrLK(img1, img2, left_points, right_points, vstatus, verror);

	std::vector<Point2f> right_points_to_find;
	std::vector<int> right_points_to_find_back_index;

	for (unsigned int i = 0; i < vstatus.size(); i++)
	{
		if (vstatus[i] && verror[i] < 8.0)
		{
			right_points_to_find_back_index.push_back(i);
			right_points_to_find.push_back(right_points[i]);
		}
		else {
			vstatus[i] = 0;
		}
	}

	Mat right_points_to_find_flat = Mat(right_points_to_find).reshape(1, right_points_to_find.size());

	std::vector<Point2f> right_features;
	KeypointsToPoints(right_keypoints, right_features);

	Mat right_features_flat = Mat(right_features).reshape(1, right_features.size());

	BFMatcher matcher(CV_L2);
	std::vector<std::vector<DMatch> >nearest_neighbors;
	matcher.radiusMatch(right_points_to_find_flat, right_features_flat, nearest_neighbors, 2.0f);

	std::set<int> found_in_right_points;
	for (int i = 0; i < nearest_neighbors.size(); i++)
	{
		DMatch _m;
		if (nearest_neighbors[i].size() == 1)
		{
			_m = nearest_neighbors[i][0];
		}
		else if (nearest_neighbors[i].size() > 1)
		{
			double ratio = nearest_neighbors[i][0].distance / nearest_neighbors[i][1].distance;

			if (ratio < 0.3) {
				_m = nearest_neighbors[i][0];
			}
			else {
				continue;
			}
		}
		else
		{
			continue;
		}

		if (found_in_right_points.find(_m.trainIdx) == found_in_right_points.end())
		{
			_m.queryIdx = right_points_to_find_back_index[_m.queryIdx];
			matches.push_back(_m);
			found_in_right_points.insert(_m.trainIdx);
		}
	}
}


void CalMatchesThroughMatcher(const Mat &img1, const Mat &img2,
	std::vector<KeyPoint> &left_keypoints, std::vector<KeyPoint> &right_keypoints,
	std::vector<DMatch> &matches)
{
	SurfDescriptorExtractor surfDesc;
	Mat descriptros1, descriptros2;
	surfDesc.compute(img1, left_keypoints, descriptros1);
	surfDesc.compute(img2, right_keypoints, descriptros2);

	//BruteForceMatcher<L2<float>> testMatcher;
	BFMatcher testMatcher(NORM_L2, true);//2.4.10, calling parameters might be different 

	testMatcher.match(descriptros1, descriptros2, matches);
	if (matches.size() > 9) {
		std::nth_element(matches.begin(), matches.begin() + 8, matches.end());//find 9 matches with smallest distances 
		matches.erase(matches.begin() + 9, matches.end());//erase the less larger distance
	}
	
}

void GetBestMatches(const std::vector<DMatch> &matches1, const std::vector<DMatch> &matches2, std::vector<DMatch> &bestMatch)
{
	for (int i = 0; i < matches1.size(); i++)
	{
		for (int j = 0; j < matches2.size(); j++)
		{
			if (matches1[i].trainIdx == matches2[j].trainIdx && matches1[i].queryIdx == matches2[j].queryIdx)
			{
				bestMatch.push_back(matches1[i]);
				break;
			}
		}
	}
}

int tencentrun()
{
	Mat img1_orin = imread(".\\topS\\Segment Fault8\\0.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	Mat img2_orin = imread(".\\topS\\Segment Fault8\\1.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	Mat img1 = img1_orin;
	Mat img2 = img2_orin;

	resize(img1_orin, img1, Size(640, 360), 0, 0, CV_INTER_LINEAR);
	resize(img2_orin, img2, Size(640, 360), 0, 0, CV_INTER_LINEAR);

	int minHessian = 60;

	SurfFeatureDetector detector(minHessian);

	//SiftFeatureDetector detector(minHessian);

	std::vector<KeyPoint> left_keypoints, right_keypoints;
	std::vector<DMatch> matches;
	detector.detect(img1, left_keypoints);
	detector.detect(img2, right_keypoints);

	std::vector<DMatch> matches1, matches2;
	//CalMatchesThroughOpticalFlk(img1, img2, left_keypoints, right_keypoints, matches);

	if (left_keypoints.size() < 10 || right_keypoints.size() < 10) {
		std::cout << "Not enough key points." << std::endl;
		return -1;
	}

	CalMatchesThroughMatcher(img1, img2, left_keypoints, right_keypoints, matches);

	//GetBestMatches(matches1, matches2, matches);
	std::vector<Point2f> imgpts1, imgpts2;
	for (unsigned int i = 0; i < matches.size(); i++)
	{
		imgpts1.push_back(left_keypoints[matches[i].queryIdx].pt);
		imgpts2.push_back(right_keypoints[matches[i].trainIdx].pt);
	}

	std::vector<DMatch> bestMatches;

	if (matches.size() > 7)
	{
		std::vector<uchar> status;
		Mat F = findFundamentalMat(imgpts1, imgpts2, FM_RANSAC, 0.01f, 0.9f, status);

		std::cout << F << std::endl;

		for (int i = 0; i < matches.size(); i++)
		{
			if (status[i])
			{
				bestMatches.push_back(matches[i]);
			}
		}
		
		//float value[] = { imgpts1[0].x, imgpts2[0].y };
		//Mat t(3, 1, CV_32FC1, value);
		//std::cout << t;
		//Mat R = F * t;


		//std::vector<unsigned char> inlinerMask(imgpts1.size());
		//Mat H = findHomography(imgpts1, imgpts2, CV_RANSAC, 1.0f, inlinerMask);
		//std::cout << H << std::endl;
		
	}
	else {
		bestMatches = matches;
	}
	std::sort(bestMatches.begin(), bestMatches.end());

	if (max(bestMatches[0].distance, bestMatches[1].distance) >= maxError)
	{
		std::cout << "Can't get a goodResult!" << std::endl;

	}
	else {
		std::cout << "Get a pretty good result!" << std::endl;
	}
	float angle, dis;

	tCalRotationAndTranslation(left_keypoints, right_keypoints, bestMatches, angle, dis);
	printf("rotation angle is %f\n", angle);
	//printf("translation is %f\n", dis);

	
	//std::vector<Point2f> newImgPt1(bestMatches.size());
	//std::vector<Point2f> newImgPt2(bestMatches.size());

	//for (int i = 0; i < bestMatches.size();i++)
	//{
	//newImgPt1[i] = left_keypoints[bestMatches[i].queryIdx].pt;
	//newImgPt2[i] = right_keypoints[bestMatches[i].trainIdx].pt;
	//}
	//DrawEpiLines(img1, img2, newImgPt1, newImgPt2);
	

	Mat img_matches;
	drawMatches(img1, left_keypoints, img2, right_keypoints, bestMatches, img_matches, Scalar::all(-1), Scalar::all(-1),
		vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	//DrawEpiLines(img1, img2, imgpts1, imgpts2);
	imshow("Good Matches & Object detection", img_matches);

	waitKey(0);
	return 0;
}
*/
///////
///////
void t1();
void t2();
void t5();
void t6();
int main()
{
	//t5();
	//t6();
	//t1();
	//t2();
	//tencentrun();
	Val2();
	system("pause");
    return 0;
}

void t1() {
	//Mat img1_orin = imread(".\\topS\\Segment Fault8\\0.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	//Mat img2_orin = imread(".\\topS\\Segment Fault8\\1.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	Mat img1_orin = imread(".\\1\\2_65107_o.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	Mat img2_orin = imread(".\\1\\2_65107_a.jpg", CV_LOAD_IMAGE_GRAYSCALE);
	EstRXY e(img1_orin.cols, img1_orin.rows);
	double r;
	Point2f xy_v;
	TickMeter time;
	time.start();
	int flag = e.estRXYVector(img1_orin, img2_orin, r, xy_v);
	time.stop();
	printf("Time consuming: %fs", time.getTimeSec());
	if (flag>=0) {
		std::printf("\nFlag:%d\nrotation:%lf\n", flag, r);
		std::printf("(x,y): (%f, %f)\n", xy_v.x, xy_v.y);
	}else{
		printf("\nFalse detection.\n");
	}
}
void t2() {
	ImageSet imgs;
	//groupFiles(".\\topS\\Segment Fault8", imgs);
	groupFiles(".\\1", imgs);
	readDirectorySeg(imgs);
	int count;
	count = 0;
	for (size_t i = 0; i < imgs.fc.size()-1; i++)
	{
		EstRXY e(imgs.fc[0].cols, imgs.fc[0].rows);
		double r;
		Point2f xy_v;
		TickMeter time;
		time.start();
		int flag = e.estRXYVector(imgs.fc[i], imgs.fc[i+1], r, xy_v);
		i += 1;
		time.stop();
		count += 1;
		printf("\n%d Time consuming: %fs", count, time.getTimeSec());
		if (flag >= 0) {
			std::printf("\nFlag:%d\nrotation:%lf\n", flag, r);
			std::printf("(x,y): (%f, %f)\n", xy_v.x, xy_v.y);
		}
		else {
			printf("\nFalse detection.\n");
		}
	}
}
void t5() {
	TickMeter t1;
	ImageSet imgs1, imgs2;
	groupFiles(".\\图片样本\\第一组\\第一点", imgs1);
	groupFiles(".\\图片样本\\第一组\\第二点", imgs2);
	readDirectorySeg(imgs1);
	readDirectorySeg(imgs2);

	SYSTEMTIME st;
	char sdir[1024], n1[1024], n2[1024];
	GetLocalTime(&st);
	sprintf(sdir, ".\\图片样本\\第一组\\TC_%4d%02d%02d_%02d%02d summary",
		st.wYear, st.wMonth, st.wDay,
		st.wHour, st.wMinute);
	FILE* f = fopen(sdir, "w");
	int size1 = imgs1.fc.size();
	int size2 = imgs2.fc.size();
	Mat am, ima, amh;
	double r, x, y;
	int pos;
	t1.start();
	int count = 0;
	for (int i = 0;i < size1;i++) {
		for (int j = 0;j < size2;j++) {
			pos = imgs1.dir[i].rfind("\\");
			pos = imgs1.dir[i].copy(n1, imgs1.dir[i].length() - pos - 5, pos + 1);
			n1[pos] = '\0';
			pos = imgs2.dir[j].rfind("\\");
			pos = imgs2.dir[j].copy(n2, imgs2.dir[j].length() - pos - 5, pos + 1);
			n2[pos] = '\0';
			count += 1;


			/*Highlight: usage of EstRXY*/
			EstRXY e(imgs1.fc[i].cols,imgs1.fc[i].rows);
			//if(e.estRXY(imgs2.fc[j],imgs1.fc[i],r/*,x,y*/)){


			//	fprintf(f, "%s %s 14 14 %lf\n",
			//		n1,n2,
			//		r-(atoi(n2) - atoi(n1))/*,
			//		x-14,y-14*/);
			//	std::printf("%d %s %s %d %d %lf\n",
			//		count,n1,n2,
			//		i, j,r/*, x, y*/);
			//}
			//else {
			//	std::printf("%d %s %s %d %d false detection\n\n",
			//		count,n1,n2,i, j);
			//}
		}
	}
	t1.stop();
	std::printf("\nTime: %lf\ns/p: %lf\n", t1.getTimeSec(), t1.getTimeSec() / (size1*size2));
	std::fclose(f);
}

void t6() {
	TickMeter t1;
	ImageSet imgs1, imgs2;
	groupFiles(".\\图片样本\\第二组\\第一点", imgs1);
	groupFiles(".\\图片样本\\第二组\\第二点", imgs2);
	readDirectorySeg(imgs1);
	readDirectorySeg(imgs2);

	SYSTEMTIME st;
	char sdir[1024], n1[1024], n2[1024];
	GetLocalTime(&st);
	sprintf(sdir, ".\\图片样本\\第二组\\TC_%4d%02d%02d_%02d%02d summary",
		st.wYear, st.wMonth, st.wDay,
		st.wHour, st.wMinute);
	FILE* f = fopen(sdir, "w");
	int size1 = imgs1.fc.size();
	int size2 = imgs2.fc.size();
	double r, x, y;
	int pos;
	t1.start();
	int count = 0;
	for (int i = 0;i < size1;i++) {
		for (int j = 0;j < size2;j++) {
			pos = imgs1.dir[i].rfind("\\");
			pos = imgs1.dir[i].copy(n1, imgs1.dir[i].length() - pos - 5, pos + 1);
			n1[pos] = '\0';
			pos = imgs2.dir[j].rfind("\\");
			pos = imgs2.dir[j].copy(n2, imgs2.dir[j].length() - pos - 5, pos + 1);
			n2[pos] = '\0';
			count += 1;


			/*Highlight: usage of EstRXY*/
			EstRXY e(imgs1.fc[i].cols, imgs1.fc[i].rows);
			//if (e.estRXY(imgs2.fc[j], imgs1.fc[i], r/*, x, y*/)) {


			//	fprintf(f, "%s %s 21 7 %lf\n",
			//		n1,n2,r - (atoi(n2) - atoi(n1))/*,
			//		x - 21, y - 7*/);
			//	std::printf("%d %s %s %d %d %lf\n",
			//		count,n1,n2,i, j,r/*, x, y*/);
			//}
			//else {
			//	std::printf("%d %s %s %d %d false detection\n\n",
			//		count,n1,n2,i, j);
			//}
		}
	}
	t1.stop();
	std::printf("\nTime: %lf\ns/p: %lf\n", t1.getTimeSec(), t1.getTimeSec() / (size1*size2));
	std::fclose(f);
}

