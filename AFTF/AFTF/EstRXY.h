#pragma once
#ifndef ESTRXY_H
#define ESTRXY_H

#include "stdafx.h"

class EstRXY {
public:

	/*default use a boundary metric in validation with
	ransacT=4
	boundary value for rotation  20бу*/
	EstRXY(int img_cols, int img_rows);

	/*ransacThresh: allow more points used in estimation with a bigger value, usually range from 1 - 9
	r_bound: absolute boundary for the rotation estimation*/
	EstRXY(int img_cols, int img_rows, double ransacThresh,
		double r_bound = (20.0));

	/*Must be called after estRXYVector function*/
	double getRotationEst();

	/*Must be called after estRXYVector function*/
	cv::Point2f getXYPixelEstWithinSamePlain();

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
	int estRXYVector(cv::Mat& src_in, cv::Mat& dst_in,
		double& e_rotation_out, Point2f& xy_vector_out);
	
private:
	int RegionMaximum;//define the maximum region matched point count
	int MinimunKeyPointsRequired;//define the minimum key point number required
	int ReliableMatchedPointsNumber;//define the matched point number that could get a reliable estimation
	double MaximunDistanceVariance;//define the maximun distance variance
	double ReliableDistributionVariance;//define the distribution variance that could get a reliable estimation
	
	double r;//rotation difference estimation from src to dst
	cv::Point2f center;//the center coordinate of the images
	cv::Point2f xy_vector_within_the_same_plain;//(x,y) vector of pixel difference between two input images
	double rbound;//rotation boundary value to validate the final estimation
	double ransacT;//ransac threshold
	
	/*convert img1 and img2 to gray images and save the data back in original matrix,
	if img1 and img2 are  BGR images. And perform equalHist on img1 and img2 if they are gray image.*/
	void normImages(cv::Mat& img1, cv::Mat& img2);

	/*use SURF descriptor to find loosely matched points
	Return values:
	-2: failed to find transform matrix
    -1: estimation out of boundary
	0: the estimation is not reliable
	1: the estimation might be not reliable
	2: the estimation is reliable*/
	int useSurfDescriptors(cv::Mat& src_in, cv::Mat& dst_in);

	/*adjust cf1 point set*/
	void adjustMatchedPoints(vector<Point2f>& cf1);

	/*filt out the ummatched points*/
	void  filtOutInvalidPoints(std::vector<Point2f>& p1in, std::vector<Point2f>& p2in, 
		std::vector<Point2f>& p1out, std::vector<Point2f>& p2out, 
		std::vector<unsigned char>& match_mask);
	
	/*return true if the transformMatrix is a 2*3 matrix and calculate the estimations*/
	int calEstRXY(vector<Point2f>&cf1, vector<Point2f>&cf2);

	/*return true if and only if the rotation estimation is within the boundary*/
	bool validateEstimation();

};

#endif // !ESTRXY_H