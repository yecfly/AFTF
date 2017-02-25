#include "stdafx.h"
#include "EstRXY.h"

using namespace std;
using namespace cv;


EstRXY::EstRXY(int img_cols, int img_rows) {
	this->RegionMaximum = 8;
	this->MinimunKeyPointsRequired = 300;
	this->MaximunDistanceVariance = 100;
	this->ReliableMatchedPointsNumber = 1000;
	this->ReliableDistributionVariance = 7.93;
	this->ransacT = 4;
	this->rbound = 20;
	this->center.x = (float)img_cols / 2;
	this->center.y = (float)img_rows / 2;
}

EstRXY::EstRXY(int img_cols, int img_rows, double ransacThresh, double r_bound/*, double x_bound, double y_bound*/) {
	this->RegionMaximum = 8;
	this->MinimunKeyPointsRequired = 300;
	this->MaximunDistanceVariance = 100; 
	this->ReliableMatchedPointsNumber = 1000;
	this->ReliableDistributionVariance = 7.93;
	this->ransacT = ransacThresh;
	this->rbound = r_bound;
	this->center.x = (float)img_cols / 2;
	this->center.y = (float)img_rows / 2;
}

double EstRXY::getRotationEst() {
	return this->r;
}

cv::Point2f EstRXY::getXYPixelEstWithinSamePlain() {
	return this->xy_vector_within_the_same_plain;
}

void EstRXY::normImages(Mat& img1, Mat& img2) {
	if (img1.channels() == 3) {
		cvtColor(img1, img1, CV_BGR2GRAY);
	}
	if (img2.channels() == 3) {
		cvtColor(img2, img2, CV_BGR2GRAY);
	}
	if (img1.channels() == 1) {
		equalizeHist(img1, img1);
	}
	if (img2.channels() == 1) {
		equalizeHist(img2, img2);
	}
}

void EstRXY::filtOutInvalidPoints(vector<Point2f>& k1in, vector<Point2f>& k2in, vector<Point2f>& k1out, vector<Point2f>& k2out, vector<unsigned char>& match_mask) {
	int size = match_mask.size();
	k1out.clear();
	k2out.clear();
	for (int i = 0;i < size;i++) {
		if (match_mask[i] > 0) {
			k1out.push_back(k1in[i]);
			k2out.push_back(k2in[i]);
		}
	}
}

void EstRXY::adjustMatchedPoints(vector<Point2f>& cf1) {
	double rot = 3.1415926535898*this->r / 180.0;
	double sinr = sin(rot);
	double cosr = cos(rot);
	Point2f tem;
	for (int i = 0;i < cf1.size();i++) {
		cf1[i] = cf1[i] - this->center;
		tem.x = cf1[i].x*cosr - sinr*cf1[i].y;
		tem.y = cf1[i].x*sinr + cf1[i].y*cosr;
		cf1[i] = tem + this->center;
	}
}

/*Estimate the rotaion, vector (x,y) between two images(should be in the same plain).
Return true with output in parameters if the estimation assume to be reasonable.
Otherwise, return false with no parameters being set.
src_in: image after movement, should be in gray scale
dst_in: image before movement, should be in gray scale
e_rotation_out: rotation estimation
xy_vector: (x,y) pixel vector estimation
Return values:
-2: failed to find transform matrix
-1: estimation out of boundary
0: the estimation is not reliable
1: the estimation might be not reliable
2: the estimation is reliable*/
int EstRXY::estRXYVector(cv::Mat& src_in, cv::Mat& dst_in,
	double& e_rotation_out, Point2f& xy_vector_out) {
	Mat am;
	normImages(src_in, dst_in);
	int flag = useSurfDescriptors(src_in, dst_in);
	if (flag>=0) {
		e_rotation_out = this->r;
		xy_vector_out = this->xy_vector_within_the_same_plain;
	}
	return flag;
}


int EstRXY::useSurfDescriptors(Mat& src_in, Mat& dst_in) {
	//int minHessian = 400;
	SurfFeatureDetector det(400);
	vector<KeyPoint> kps1, kps2;
	Mat desc1, desc2;
	det.detect(src_in, kps1);
	det.detect(dst_in, kps2);
	if (kps1.size() < this->MinimunKeyPointsRequired || kps2.size() < this->MinimunKeyPointsRequired) {
		printf("\n******************\nError: The input images are not meeting the minimun \nrequirement of key point numbers.\n");
		return -1;
	}
	SurfDescriptorExtractor ext;
	ext.compute(src_in, kps1, desc1);
	ext.compute(dst_in, kps2, desc2);

	//FlannBasedMatcher matcher;
	BFMatcher matcher(NORM_L2, true);
	vector<DMatch> matches;
	matcher.match(desc1, desc2, matches);
	int sizem = matches.size();
	vector<Point2f> cp1, cp2, cf1, cf2;
	for (int i = 0;i < sizem;i++) {
		cp1.push_back(kps1[matches[i].queryIdx].pt);
		cp2.push_back(kps2[matches[i].trainIdx].pt);
	}
	vector<unsigned char> match_mask;
	findHomography(cp1, cp2, CV_RANSAC, this->ransacT, match_mask);
	filtOutInvalidPoints(cp1, cp2, cf1, cf2, match_mask);
	
	return calEstRXY(cf1, cf2);
}

int EstRXY::calEstRXY(vector<Point2f>&cf1, vector<Point2f>&cf2) {
	Mat am;
	estimateRigidTransform(cf1, cf2, false).copyTo(am);
	if (am.cols == 3 && am.rows == 2) {
		/*calculate rotation*/
		double xcos = am.at<double>(0, 0);
		double xsin = am.at<double>(0, 1);
		double sqr = sqrt(xcos*xcos + xsin*xsin);
		double sinv = xsin / sqr;
		double cosv = xcos / sqr;
		this->r = 11215.0 - 250.05*sinv - 11215.0*cosv - 5650.0*sinv*sinv + 300.0*sinv*cosv;
		
		/*adjust cf1 for a mean vector*/
		adjustMatchedPoints(cf1);
		Scalar tmdis = sum(cf1) - sum(cf2);
		Point2f xyv = Point2f(tmdis[0] / (double)cf1.size(), tmdis[1] / (double)cf1.size());
		double x, y;

		/*calculate y axis movement*/
		x = sinv*xyv.y;
		y = cosv*xyv.x;
		this->xy_vector_within_the_same_plain.y = 0.02715 - 1.25*x - 0.01937*y
			- 0.02628*x*x + 0.004281*x*y - 0.000168*y*y
			+ 0.001796*x*x*x + (9.265e-05)*x*x*y + (1.202e-05)*x*y*y;

		/*calculate x axis movement*/
		x = sinv*xyv.x;
		y = cosv*xyv.y;
		this->xy_vector_within_the_same_plain.x = -0.04476 - 0.04489*x + 0.05912*y
			- 0.00397*x*x + 0.002685*x*y - (1.138e-05)*y*y
			+ (4.457e-05)*x*x*x - (5.631e-05)*x*x*y - (2.755e-05)*x*y*y;

		/*validate the estimation and its reliability*/
		if (validateEstimation()) {
			/*matched points analyses*/
			int size = cf1.size();
			if (size >= this->ReliableMatchedPointsNumber) {
				return 2;//reliable estimation
			}

			/*distance analyses*/
			vector<double> evad;
			evad.resize(size);
			Point2f tmd;
			for (int i = 0;i < size;i++) {
				tmd = cf1[i] - cf2[i];
				evad[i] = sqrt(tmd.dot(tmd));
			}
			Scalar tm = sum(evad);
			double mean = tm[0] / (double)size;
			double variance = 0;
			double tvar;
			for (int i = 0;i < size;i++) {
				tvar = evad[i] - mean;
				variance += sqrt(tvar*tvar);
			}
			variance = variance / (double)size;
			if (variance > this->MaximunDistanceVariance) {
				return 0;//The estimation is not reliable
			}

			/*distribution analyses*/
			vector<float> eva_distribution;//evaluate the distribution of matched points
			eva_distribution.resize(4);
			eva_distribution[0] = 0;
			eva_distribution[1] = 0;
			eva_distribution[2] = 0;
			eva_distribution[3] = 0;
			Point2f t;
			for (int i = 0; i < size; i++) {
				t = cf1[i] - this->center;
				if(t.x<=0){
					if (t.y <= 0) {
						/*upper region on the left, region 1*/
						eva_distribution[0] += 1;
					}
					else {
						/*upper region on the right, region 2*/
						eva_distribution[1] += 1;
					}
				}else {
					if (t.y <= 0) {
						/*lower region on the left, region 3*/
						eva_distribution[2] += 1;
					}
					else {
						/*lower region on the right, region 4*/
						eva_distribution[3] += 1;
					}
				}
			}
			for (int i = 0;i < 4;i++) {
				if (eva_distribution[i] > this->RegionMaximum) {
					eva_distribution[i] = this->RegionMaximum;
				}
			}
			tm = sum(eva_distribution);
			tm[1] = tm[0] / (double)size;
			tm[2] = 0;
			for (int i = 0;i < 4;i++) {
				tvar = eva_distribution[i] - tm[1];
				tm[2] = tm[2] + sqrt(tvar*tvar);
			}
			tm[2] = tm[2] / 4.0;
			if (tm[2] > this->ReliableDistributionVariance) {
				return 2;//reliable estimation
			}
			
			return 1;//the estimation might be not reliable
		}
		else {
			return -1;//estimation out of boundary
		}
	}
	return -2;//failed to find transform matrix
}

bool EstRXY::validateEstimation() {
	if (abs(this->r) < this->rbound 
		&& abs(this->xy_vector_within_the_same_plain.x) < 2*this->center.x 
		&& abs(this->xy_vector_within_the_same_plain.y) < 2*this->center.y) {
		return true;
	}
	std::printf("\nWARNING: estimation exceeds the boundary.\n");
	return false;
}
