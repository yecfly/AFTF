
#include "stdafx.h"
#include <direct.h>

using namespace std;
using namespace cv;

struct estRXY {
	double rotation;//in form of 180°, 2 means 2°
	double x;//with a unit of mm, 2 means 2mm
	double y;//with a unit of mm, 2 means 2mm
};

struct PSE {
	double mean;//mean distance between two point set
	double variance;//variance between two point set
};
//group the files in the directory to prevent over loading, return the size
//of vfs
int groupFiles(const string& directoryName, ImageSet& vfs) {
	vector<string> dir;
	dir.clear();

	struct _finddata_t s_file;
	string str = directoryName + "\\*.*";

	intptr_t h_file = _findfirst(str.c_str(), &s_file);
	if (h_file != static_cast<intptr_t>(-1.0))
	{
		//_findnext(h_file, &s_file);//pass /.
		//_findnext(h_file, &s_file);//pass /..
		while (_findnext(h_file, &s_file) == 0) {
			if (s_file.attrib == _A_ARCH) {
				dir.push_back(directoryName + "\\" + s_file.name);
			}
		}
	}
	_findclose(h_file);

	//only remain the .jpg extention files' dir
	string ext = ".jpg";
	string ext1 = ".png";
	string ext2 = ".bmp";
	vector<string>::iterator it;
	it = dir.begin();
	while (it != dir.end()) {
		string y = *it;
		if ((y.rfind(ext) != y.length() - 4) && (y.rfind(ext1) != y.length() - 4) && (y.rfind(ext2) != y.length() - 4)) {
			it = dir.erase(it);
		}
		else {
			it++;
		}
	}

	//compare with the maximageload
	if (dir.size() > 0) {
		sort(dir.begin(), dir.end());
		vector<string>::iterator ite;
		ite = dir.begin();
		vfs.dir.assign(ite, dir.end());
	}
	return vfs.dir.size();
}

//return the image's absolute readable path to ImageSet.dir and the
//image data to ImageSet.frame
void readDirectorySeg(ImageSet& frames) {
	int size_f = frames.dir.size();
	for (int i = 0;i < size_f;i++) {
		frames.fc.push_back(imread(frames.dir[i], IMREAD_GRAYSCALE));
		equalizeHist(frames.fc[i], frames.fc[i]);
	}
}

bool calEstR(double& r, cv::Mat& transformMatrix) {
	if (transformMatrix.cols == 3 && transformMatrix.rows == 2) {
		double xcos = transformMatrix.at<double>(0, 0);
		double xsin = transformMatrix.at<double>(0, 1);
		double sqr = sqrt(xcos*xcos + xsin*xsin);
		double sinv = xsin / sqr;
		double cosv = xcos / sqr;

		//e.rotation = -73.43 + 96.5*sinv + 71.9*cosv;
		//r = 22430.0 - 500.1*sinv - 22430.0*cosv - 11300.0*sinv*sinv + 600.0*sinv*cosv;//filter out 6 points
		//r = r / 2;// bias adjustment
		r = 11215.0 - 250.05*sinv - 11215.0*cosv - 5650.0*sinv*sinv + 300.0*sinv*cosv;																					//e.rotation = 27780.0 - 518.5*sinv - 27780.0*cosv - 14000.0*sinv*sinv + 617.9*sinv*cosv;//filter out 7 points, poly21
		
		return true;
	}
	return false;
}

int calXYV(vector<Point2f>&cf1, vector<Point2f>&cf2, Point2f& center, Point2f& xyv) {
	int size = cf1.size();
	if (size == cf2.size() && size > 3) {

		float lu, ru, ll, rl, dis;
		lu = ru = ll = rl = 0;
		vector<float> eva_distribution;//evaluate the distribution of matched points
		eva_distribution.resize(4);
		eva_distribution[0] = 0;
		eva_distribution[1] = 0;
		eva_distribution[2] = 0;
		eva_distribution[3] = 0;
		Point2f t;//c1 c2 c3 c4 represent the matched points with longest distance to the center inside each region
		vector<Point2f> c, cp;
		c.resize(4);
		cp.resize(4);
		for (int i = 0; i < size; i++) {
			t = cf1[i] - center;
			if (t.x == 0 && t.y == 0) {//found the corespongding center point in img2
				xyv = center - cf2[i];
				return 0;
			}
			else if (t.x <= 0) {
				dis = t.dot(t);
				if (t.y <= 0) {
					/*upper region on the left, region 1*/
					if (dis > lu) {
						lu = dis;
						c[0] = cf1[i];
						cp[0] = cf2[i];
						eva_distribution[0] += 1;
					}
				}
				else {
					/*upper region on the right, region 2*/
					if (dis > ru) {
						ru = dis;
						c[1] = cf1[i];
						cp[1] = cf2[i];
						eva_distribution[1] += 1;
					}
				}
			}
			else {
				dis = t.dot(t);
				if (t.y <= 0) {
					/*lower region on the left, region 3*/
					if (dis > ll) {
						ll = dis;
						c[2] = cf1[i];
						cp[2] = cf2[i];
						eva_distribution[2] += 1;
					}
				}
				else {
					/*lower region on the right, region 4*/
					if (dis > rl) {
						rl = dis;
						c[3] = cf1[i];
						cp[3] = cf2[i];
						eva_distribution[3] += 1;
					}
				}
			}
		}
		vector<Point2f> v1, v2;
		for (int j = 0; j < 4; j++) {//normalize the distribution
			if (eva_distribution[j] > 3) {
				eva_distribution[j] = 3;
				v1.push_back(c[j]);
				v2.push_back(cp[j]);
			}
			else if (eva_distribution[j]>0) {
				v1.push_back(c[j]);
				v2.push_back(cp[j]);
			}
		}
		int psize = v1.size();
		vector<Point2f> xyr;
		float xr, yr, t12, t13, xm2, ym2, xm3, ym3, xv21, yv21;
		if (psize > 2) {
			for (int s1 = 0; s1 < psize; s1++) {
				for (int s2 = s1 + 1; s2 < psize; s2++) {
					for (int s3 = s2 + 1; s3 < psize; s3++) {
						xm2 = v2[s2].x - v2[s1].x;
						xm3 = v2[s3].x - v2[s1].x;
						ym2 = v2[s2].y - v2[s1].y;
						ym3 = v2[s3].y - v2[s1].y;
						xv21 = v2[s1].x;
						yv21 = v2[s1].y;
						t12 = (center - v1[s1]).dot(v1[s2] - v1[s1]);
						t13 = (center - v1[s1]).dot(v1[s3] - v1[s1]);
						if (ym2 == 0 && ym3 != 0 && xm2 != 0) {
							xr = t12 / xm2 + xv21;
							yr = (t13 - xm3*(xr - xv21)) / ym3 + yv21;
							xyr.push_back(Point2f(xr, yr));
						}
						else if (ym3 == 0 && ym2 != 0 && xm3 != 0) {
							xr = t13 / xm3 + xv21;
							yr = (t12 - xm2*(xr - xv21)) / ym2 + yv21;
							xyr.push_back(Point2f(xr, yr));
						}
						else if (ym2 != 0 && (ym2*xm3 - ym3*xm2) != 0) {
							xr = ym2*(t13 + xm3*xv21 + ym3*yv21) - ym3*(t12 + xm2*xv21 + ym2*yv21);
							xr = xr / (ym2*xm3 - ym3*xm2);
							yr = (t12 + xm2*xv21 + ym2*yv21 - xm2*xr) / ym2;
							xyr.push_back(Point2f(xr, yr));
						}
					}
				}
			}
			Scalar teva = sum(xyr);
			if (xyr.size() > 0) {
				xyv.x = (float)teva[0] / (float)xyr.size();
				xyv.y = (float)teva[1] / (float)xyr.size();
				xyv = center - xyv;
			}
			else {
				return -1;//unreliable case
			}
			teva = sum(eva_distribution);
			if (teva[0] > 3 * 3) {
				return 2;//reliable case
			}
			else if (teva[0] > 2 * 3) {
				return 1;//less reliable case
			}
			else {
				return 0;//might be not reliable
			}
		}
		else {
			return -1;//unreliable case
		}
	}
	else {
		return -1;//unreliable case
	}
}

void adjustV1(vector<Point2f>& cf1, Point2f& center, double& r) {
	double pi = 3.1415926535898;
	double rot = 3.1415926535898*r / 180.0;
	double sinr = sin(rot);
	double cosr = cos(rot);
	Point2f tem;
	for (int i = 0;i < cf1.size();i++) {
		cf1[i] = cf1[i] - center;
		tem.x = cf1[i].x*cosr - sinr*cf1[i].y;
		tem.y = cf1[i].x*sinr + cf1[i].y*cosr;
		cf1[i] = tem + center;
	}
}

double d2p(Point2f& p1, Point2f& p2) {
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

void PointSetEva(vector<Point2f>&set1, vector<Point2f>&set2, PSE& pse) {
	vector<double> evad;
	int size = set1.size();
	evad.resize(size);
	for (int i = 0;i < size;i++) {
		evad[i] = d2p(set1[i], set2[i]);
	}
	Scalar tm = sum(evad);
	pse.mean = tm[0] / size;
	pse.variance = 0;
	double t;
	for (int i = 0;i < size;i++) {
		t = evad[i] - pse.mean;
		pse.variance += sqrt(t*t);
	}
	pse.variance = pse.variance / size;
}

void filterGM(vector<DMatch>& in, vector<DMatch>& out, vector<unsigned char>& match_mask) {
	int size = match_mask.size();
	out.clear();
	for (int i = 0;i < size;i++) {
		if (match_mask[i] > 0) {
			out.push_back(in[i]);
		}
	}
}

void filterP(vector<Point2f>& k1in, vector<Point2f>& k2in, vector<Point2f>& k1out, vector<Point2f>& k2out, vector<unsigned char>& match_mask) {
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

void useSurfDes(Mat& src, Mat& dst, vector<KeyPoint>& kps1, vector<KeyPoint>& kps2,
	vector<Point2f>& ckp1, vector<Point2f>& ckp2, vector<DMatch>& gm) {
	SurfFeatureDetector det(400);
	//vector<KeyPoint> kps1, kps2;
	Mat desc1, desc2, ima;
	det.detect(src, kps1);
	det.detect(dst, kps2);
	SurfDescriptorExtractor ext;
	ext.compute(src, kps1, desc1);
	ext.compute(dst, kps2, desc2);

	//FlannBasedMatcher matcher;
	BFMatcher matcher(NORM_L2, true);
	vector<DMatch> matches;
	matcher.match(desc1, desc2, matches);
	ckp1.clear();
	ckp2.clear();
	int sizem = matches.size();
	for (int i = 0;i < sizem;i++) {
		ckp1.push_back(kps1[matches[i].queryIdx].pt);
		ckp2.push_back(kps2[matches[i].trainIdx].pt);
		gm.push_back(matches[i]);
	}
	/*drawMatches(src, kps1, dst, kps2, gm, ima, Scalar::all(-1), Scalar::all(-1),
	vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
	resize(ima, ima, Size(1920, 540));
	imshow("Good Matches", ima);
	waitKey();*/
	std::printf("\n%d/%d\n", ckp1.size(), kps1.size());
}

int calEstRXY(vector<Point2f>&cf1, vector<Point2f>&cf2,
	double& r, Point2f& xyr, Point2f& center, double& disvar) {
	Mat am;
	estimateRigidTransform(cf1, cf2, false).copyTo(am);
	if (am.cols == 3 && am.rows == 2) {
		/*calculate rotation*/
		double xcos = am.at<double>(0, 0);
		double xsin = am.at<double>(0, 1);
		double sqr = sqrt(xcos*xcos + xsin*xsin);
		double sinv = xsin / sqr;
		double cosv = xcos / sqr;
		r = 11215.0 - 250.05*sinv - 11215.0*cosv - 5650.0*sinv*sinv + 300.0*sinv*cosv;

		/*adjust cf1 for a mean vector*/
		adjustV1(cf1, center, r);
		Scalar tmdis = sum(cf1) - sum(cf2);
		Point2f xyv = Point2f(tmdis[0] / (double)cf1.size(), tmdis[1] / (double)cf1.size());
		double x, y;

		/*calculate y axis movement*/
		x = sinv*xyv.y;
		y = cosv*xyv.x;
		xyr.y = 0.02715 - 1.25*x - 0.01937*y
			- 0.02628*x*x + 0.004281*x*y - 0.000168*y*y
			+ 0.001796*x*x*x + (9.265e-05)*x*x*y + (1.202e-05)*x*y*y;

		/*calculate x axis movement*/
		x = sinv*xyv.x;
		y = cosv*xyv.y;
		xyr.x = -0.04476 - 0.04489*x + 0.05912*y
			- 0.00397*x*x + 0.002685*x*y - (1.138e-05)*y*y
			+ (4.457e-05)*x*x*x - (5.631e-05)*x*x*y - (2.755e-05)*x*y*y;

	    /*validate the estimation and its reliability*/
		/*matched points analyses*/
		int size = cf1.size();

		vector<float> eva_distribution;//evaluate the distribution of matched points
		eva_distribution.resize(4);
		eva_distribution[0] = 0;
		eva_distribution[1] = 0;
		eva_distribution[2] = 0;
		eva_distribution[3] = 0;
		Point2f t;
		for (int i = 0; i < size; i++) {
			t = cf1[i] - center;
			if (t.x <= 0) {
				if (t.y <= 0) {
					/*upper region on the left, region 1*/
					eva_distribution[0] += 1;
				}
				else {
					/*upper region on the right, region 2*/
					eva_distribution[1] += 1;
				}
			}
			else {
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
			if (eva_distribution[i] > 8) {
				eva_distribution[i] = 8;
			}
		}
		Scalar tm = sum(eva_distribution);
		tm[1] = tm[0] / (double)size;
		double tvar;
		disvar = 0;
		for (int i = 0;i < 4;i++) {
			tvar = eva_distribution[i] - tm[1];
			disvar = disvar + sqrt(tvar*tvar);
		}
		disvar = disvar / 4.0;
		//if (disvar < 0) {
		//	return -20;//reliable estimation
		//}

		//if (size >= 1000) {
		//	return 2;//reliable estimation
		//}

		/*distance analyses*/
		vector<double> evad;
		evad.resize(size);
		Point2f tmpd;
		for (int i = 0;i < size;i++) {
			tmpd = cf1[i] - cf2[i];
			evad[i] = sqrt(tmpd.dot(tmpd));
		}
		tm = sum(evad);
		double mean = tm[0] / (double)size;
		double variance = 0;
		for (int i = 0;i < size;i++) {
			tvar = evad[i] - mean;
			variance += sqrt(tvar*tvar);
		}
		variance = variance / (double)size;
		//if (variance > 100) {
		//	return 0;//The estimation is not reliable
		//}

		

		return 1;//the estimation might be not reliable
	}
	return -2;//failed to find transform matrix
}

void Val1() {
	TickMeter t1, t2;
	vector<string> dirs;
	dirs.push_back(".\\0623TestPic\\第一组");
	dirs.push_back(".\\0623TestPic\\第二组");
	dirs.push_back(".\\0623TestPic\\第三组");
	dirs.push_back(".\\0623TestPic\\第四组");

	for (size_t di = 0; di < dirs.size(); di++)
	{
		ImageSet imgs1;
		groupFiles(dirs[di], imgs1);
		readDirectorySeg(imgs1);

		SYSTEMTIME st;
		char logdir[1024], sdir[1024], n1[1024], n2[1024], x1[1024], y1[1024], x2[1024], y2[1024], tn[1024], name[1024];
		GetLocalTime(&st);
		sprintf(name, ".\\0623TestPic\\G%d_%4d%02d%02d_%02d%02d_Ransac%d",
			di+1,
			st.wYear, st.wMonth, st.wDay,
			st.wHour, st.wMinute, 4);
		_mkdir(name);
		sprintf(sdir, "%s\\_summary", name);
		sprintf(logdir, "%s\\_log", name);
		FILE* f = fopen(sdir, "w");
		FILE* log = fopen(logdir, "w");
		fprintf(log, "col1:计数\tcol2:图片1角度\tcol3图片2角度\tcol9:角度估值误差\tcol13:过滤后特征点集\tcol14:过滤前特征点集\tcol22:图一索引（从零开始）\tcol23(最后一列):图二索引（从零开始）\n");
		int size1 = imgs1.fc.size();
		//int size2 = imgs2.fc.size();
		Mat am, ima, amh;
		int pos, pos2;
		double xcos, xsin, ax, ay, sqr, time;
		t1.start();
		int count = 0;
		for (int i = 0;i < size1;i++) {
			for (int j = i;j < size1;j++) {
				pos = imgs1.dir[i].rfind("z");
				pos2 = imgs1.dir[i].find(".", pos);
				pos = imgs1.dir[i].copy(n1, pos2 - pos - 1, pos + 1);
				n1[pos] = '\0';
				pos = imgs1.dir[i].rfind("x");
				pos2 = imgs1.dir[i].find(",", pos);
				pos = imgs1.dir[i].copy(x1, pos2 - pos - 1, pos + 1);
				x1[pos] = '\0';
				pos = imgs1.dir[i].rfind("y");
				pos2 = imgs1.dir[i].find(",", pos);
				pos = imgs1.dir[i].copy(y1, pos2 - pos - 1, pos + 1);
				y1[pos] = '\0';

				pos = imgs1.dir[j].rfind("z");
				pos2 = imgs1.dir[j].find(".", pos);
				pos = imgs1.dir[j].copy(n2, pos2 - pos - 1, pos + 1);
				n2[pos] = '\0';
				pos = imgs1.dir[j].rfind("x");
				pos2 = imgs1.dir[j].find(",", pos);
				pos = imgs1.dir[j].copy(x2, pos2 - pos - 1, pos + 1);
				x2[pos] = '\0';
				pos = imgs1.dir[j].rfind("y");
				pos2 = imgs1.dir[j].find(",", pos);
				pos = imgs1.dir[j].copy(y2, pos2 - pos - 1, pos + 1);
				y2[pos] = '\0';

				vector<Point2f> cp1, cp2, cf1, cf2;
				vector<KeyPoint> k1, k2;
				vector<unsigned char> match_mask;
				vector<DMatch> g, gm;

				t2.reset();
				t2.start();
				useSurfDes(imgs1.fc[j], imgs1.fc[i], k1, k2, cp1, cp2, g);
				t2.stop();
				time = t2.getTimeSec();
				printf("GetSurfDes: %lfs\n", time);
				fprintf(log, "GetSurfDes: %lfs\n", time);
				count += 1;
				if (cp1.size() > 15) {
					t2.start();
					amh = findHomography(cp1, cp2, CV_RANSAC, 4, match_mask);
					filterGM(g, gm, match_mask);
					filterP(cp1, cp2, cf1, cf2, match_mask);
					t2.stop();
					std::printf("%d/%d FilterPointSets: %lfs\n", cf1.size(), cp1.size(), t2.getTimeSec() - time);
					fprintf(log, "%d/%d FilterPointSets: %lfs\n", cf1.size(), cp1.size(), t2.getTimeSec() - time);

					time = t2.getTimeSec();
					t2.start();
					am = estimateRigidTransform(cf1, cf2, false);
					t2.stop();
					printf("EstimateRigidTransform: %lfs\n", t2.getTimeSec() - time);
					fprintf(log, "EstimateRigidTransform: %lfs\n", t2.getTimeSec() - time);
					cout << am << endl << endl;
					/*drawMatches(imgs2.fc[j], k1, imgs1.fc[i], k2, gm, ima, Scalar::all(-1), Scalar::all(-1),
						vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
					sprintf(tn, "%s\\%03d_%s_%s.jpg",
						name, count, n1, n2);
					imwrite(tn, ima);
					resize(ima, ima, Size(1920, 540));
					imshow("Good Matches", ima);
					waitKey();*/

					if (am.cols < 3 || am.rows < 2) {
						std::printf("%d %d %d %d failed to find am\n\n", di + 1,
							count, i, j);
						std::fprintf(log,"%d %d %d %d failed to find am\n\n", di + 1,
							count, i, j);
					}
					else {
						double e=-1000;
						//calEstR(e, am);
						double distribution_var = -1000;
						Point2f xyv, xyr;
						Point2f center = Point2f(imgs1.fc[i].cols / 2, imgs1.fc[i].rows / 2);
						
						time = t2.getTimeSec();
						t2.start();
						int flag_cerxy=calEstRXY(cf1, cf2, e, xyr, center,distribution_var);
						t2.stop();
						printf("RXYEstimation: %lfs\n", t2.getTimeSec() - time);
						fprintf(log, "RXYEstimation: %lfs\n", t2.getTimeSec() - time);

						time = t2.getTimeSec();
						t2.start();
						//pos = calXYV(cf1, cf2, center, xyv);
						adjustV1(cf1, center, e);
						Scalar tmd = sum(cf1) - sum(cf2);
						xyv.x = tmd[0] / (double)cf1.size();
						xyv.y = tmd[1] / (double)cf1.size();
						t2.stop();
						printf("calXYV: %lfs\n", t2.getTimeSec() - time);
						fprintf(log, "calXYV: %lfs\n", t2.getTimeSec() - time);

						time = t2.getTimeSec();
						PSE pe;
						t2.start();
						PointSetEva(cf1, cf2, pe);
						t2.stop();
						printf("PointSetsEva: %lfs\n", t2.getTimeSec() - time);
						fprintf(log, "PointSetsEva: %lfs\n", t2.getTimeSec() - time);

						if (abs(e) <= 20 && abs(xyv.x) <= 2*center.x && abs(xyv.y) <= 2*center.y) {
							xcos = am.at<double>(0, 0);
							xsin = am.at<double>(0, 1);
							ax = am.at<double>(0, 2);
							ay = am.at<double>(1, 2);
							sqr = sqrt(xcos*xcos + xsin*xsin);
							fprintf(f, "%d %s %s %s %s %s %s %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf\n",
								count, n1, n2, x1, x2, y1, y2, 
								e, e - (atof(n2) - atof(n1)),
								xyv.x, xyv.y, pos,
								cf1.size(), cp1.size(),
								pe.mean, pe.variance, ax, ay, xcos / sqr, xsin / sqr, sqr, i, j, distribution_var, xyr.x, xyr.y);
							fprintf(log, "%d %s %s %s %s %s %s %lf %lf %lf %lf %d %d %d %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf\n",
								count, n1, n2, x1, x2, y1, y2,
								e, e - (atof(n2) - atof(n1)),
								xyv.x, xyv.y, pos,
								cf1.size(), cp1.size(),
								pe.mean, pe.variance, ax, ay, xcos / sqr, xsin / sqr, sqr, i, j, distribution_var, xyr.x, xyr.y);
							std::printf("%d %d %d %d %lf\n", di + 1,
								count, i, j, e - (atof(n2) - atof(n1)));
						}
						else {
							std::printf("%d %d %d %d boundary failure\n\n", di + 1,
								count, i, j);
							std::fprintf(log, "%d %d %d %d boundary failure\n\n", di + 1,
								count, i, j);
						}
					}
				}
				else {
					std::printf("%d %d %d %d failed to find am due to insufficient keypoints\n\n",
						di + 1, count, i, j);
					std::fprintf(log, "%d %d %d %d failed to find am due to insufficient keypoints\n\n",
						di + 1, count, i, j);
				}
			}
		}
		t1.stop();
		std::printf("\nTime: %lf\ns/p: %lf\n", t1.getTimeSec(), t1.getTimeSec() / count);
		std::fprintf(log, "\nTime: %lfs \t seconds/per pair: %lf\n", t1.getTimeSec(), t1.getTimeSec() / count);
		std::fclose(f);
		std::fclose(log);
	}
}

void Val2() {
	TickMeter t1, t2;
	vector<string> dirs;
	dirs.push_back(".\\0623TestPic\\第一组");
	dirs.push_back(".\\0623TestPic\\第二组");
	dirs.push_back(".\\0623TestPic\\第三组");
	dirs.push_back(".\\0623TestPic\\第四组");

	for (size_t di = 0; di < dirs.size(); di++)
	{
		ImageSet imgs1;
		groupFiles(dirs[di], imgs1);
		readDirectorySeg(imgs1);

		SYSTEMTIME st;
		char logdir[1024], sdir[1024], n1[1024], n2[1024], x1[1024], y1[1024], x2[1024], y2[1024], tn[1024], name[1024];
		GetLocalTime(&st);
		sprintf(name, ".\\0623TestPic\\G%d_%4d%02d%02d_%02d%02d_Ransac%d",
			di + 1,
			st.wYear, st.wMonth, st.wDay,
			st.wHour, st.wMinute, 4);
		_mkdir(name);
		sprintf(sdir, "%s\\_summary", name);
		sprintf(logdir, "%s\\_log", name);
		FILE* f = fopen(sdir, "w");
		FILE* log = fopen(logdir, "w");
		fprintf(log, "col1:计数\tcol2:图片1角度\tcol3图片2角度\tcol9:角度估值误差\tcol13:过滤后特征点集\tcol14:过滤前特征点集\tcol22:图一索引（从零开始）\tcol23(最后一列):图二索引（从零开始）\n");
		int size1 = imgs1.fc.size();
		//int size2 = imgs2.fc.size();
		Mat am, ima, amh,m1,m2,m3,m4;
		int pos, pos2, c1, c2,c3,c4;
		double xcos, xsin, ax, ay, sqr, time;
		t1.start();
		int count = 0;
		for (int i = 0;i < size1;i++) {
			for (int j = i;j < size1;j++) {
				pos = imgs1.dir[i].rfind("z");
				pos2 = imgs1.dir[i].find(".", pos);
				pos = imgs1.dir[i].copy(n1, pos2 - pos - 1, pos + 1);
				n1[pos] = '\0';
				pos = imgs1.dir[i].rfind("x");
				pos2 = imgs1.dir[i].find(",", pos);
				pos = imgs1.dir[i].copy(x1, pos2 - pos - 1, pos + 1);
				x1[pos] = '\0';
				pos = imgs1.dir[i].rfind("y");
				pos2 = imgs1.dir[i].find(",", pos);
				pos = imgs1.dir[i].copy(y1, pos2 - pos - 1, pos + 1);
				y1[pos] = '\0';

				pos = imgs1.dir[j].rfind("z");
				pos2 = imgs1.dir[j].find(".", pos);
				pos = imgs1.dir[j].copy(n2, pos2 - pos - 1, pos + 1);
				n2[pos] = '\0';
				pos = imgs1.dir[j].rfind("x");
				pos2 = imgs1.dir[j].find(",", pos);
				pos = imgs1.dir[j].copy(x2, pos2 - pos - 1, pos + 1);
				x2[pos] = '\0';
				pos = imgs1.dir[j].rfind("y");
				pos2 = imgs1.dir[j].find(",", pos);
				pos = imgs1.dir[j].copy(y2, pos2 - pos - 1, pos + 1);
				y2[pos] = '\0';

				vector<Point2f> cp1, cp2, cf1, cf2;
				vector<KeyPoint> k1, k2;
				vector<unsigned char> match_mask;
				vector<DMatch> g, gm;

				t2.reset();
				t2.start();
				Mat r1(imgs1.fc[j], Rect(0, 0, 640, 360));
				r1.copyTo(m1);
				Mat b1(imgs1.fc[i], Rect(0, 0, 640, 360));
				b1.copyTo(m2);
				useSurfDes(m1,m2, k1, k2, cp1, cp2, g);
				c1 = k1.size();
				c2 = k2.size();
				c3 = cp1.size();
				c4 = cp2.size();

				Mat r2(imgs1.fc[j], Rect(640, 0, 640, 360));
				r2.copyTo(m1);
				Mat b2(imgs1.fc[i], Rect(640, 0, 640, 360));
				b2.copyTo(m2);
				useSurfDes(m1, m2, k1, k2, cp1, cp2, g);
				c1 += k1.size();
				c2 += k2.size();
				c3 += cp1.size();
				c4 += cp2.size();

				Mat r3(imgs1.fc[j], Rect(0, 360, 640, 360));
				r3.copyTo(m1);
				Mat b3(imgs1.fc[i], Rect(0, 360, 640, 360));
				b3.copyTo(m2);
				useSurfDes(m1, m2, k1, k2, cp1, cp2, g);
				c1 += k1.size();
				c2 += k2.size();
				c3 += cp1.size();
				c4 += cp2.size();

				Mat r4(imgs1.fc[j], Rect(640, 360, 640, 360));
				r4.copyTo(m1);
				Mat b4(imgs1.fc[i], Rect(640, 360, 640, 360));
				b4.copyTo(m2);
				useSurfDes(m1, m2, k1, k2, cp1, cp2, g);
				c1 += k1.size();
				c2 += k2.size();
				c3 += cp1.size();
				c4 += cp2.size();
				fprintf(log, "%d %d %d %d\n", c1, c2, c3, c4);

				useSurfDes(imgs1.fc[j], imgs1.fc[i], k1, k2, cp1, cp2, g);
				fprintf(log, "%d %d %d %d\n\n\n", k1.size(), k2.size(), cp1.size(), cp2.size());
				
				t2.stop();
			}
		}
		t1.stop();
		std::printf("\nTime: %lf\ns/p: %lf\n", t1.getTimeSec(), t1.getTimeSec() / count);
		std::fprintf(log, "\nTime: %lfs \t seconds/per pair: %lf\n", t1.getTimeSec(), t1.getTimeSec() / count);
		std::fclose(f);
		std::fclose(log);
	}
}
//void Val2() {
//	TickMeter t1, t2;
//	ImageSet imgs1, imgs2;
//	groupFiles(".\\0623TestPic\\第二组\\第一点", imgs1);
//	groupFiles(".\\0623TestPic\\第二组\\第二点", imgs2);
//	readDirectorySeg(imgs1);
//	readDirectorySeg(imgs2);
//
//	SYSTEMTIME st;
//	char sdir[1024], n1[1024], n2[1024], x1[1024], y1[1024], x2[1024], y2[1024], tn[1024], name[1024];
//	GetLocalTime(&st);
//	sprintf(name, ".\\0623TestPic\\G2_%4d%02d%02d_%02d%02d_Ransac%d",
//		st.wYear, st.wMonth, st.wDay,
//		st.wHour, st.wMinute, 4);
//	_mkdir(name);
//	sprintf(sdir, "%s\\_summary", name);
//	FILE* f = fopen(sdir, "w");
//	//fprintf(f, "图片1角度\t图片2角度\t角度估值误差\tx估值误差\ty估值误差\t过滤后特征点集\t过滤前特征点集\t角度估值\tx估值\ty估值\n");
//	int size1 = imgs1.fc.size();
//	int size2 = imgs2.fc.size();
//	Mat am, ima, amh;
//	int pos;
//	double xcos, xsin, ax, ay, sqr, time;
//	t1.start();
//	int count = 0;
//	for (int i = 0;i < size1;i++) {
//		for (int j = 0;j < size2;j++) {
//			pos = imgs1.dir[i].rfind("\\");
//			pos = imgs1.dir[i].copy(n1, imgs1.dir[i].length() - pos - 5, pos + 1);
//			n1[pos] = '\0';
//			pos = imgs2.dir[j].rfind("\\");
//			pos = imgs2.dir[j].copy(n2, imgs2.dir[j].length() - pos - 5, pos + 1);
//			n2[pos] = '\0';
//
//			vector<Point2f> cp1, cp2, cf1, cf2;
//			vector<KeyPoint> k1, k2;
//			vector<unsigned char> match_mask;
//			vector<DMatch> g, gm;
//
//			t2.reset();
//			t2.start();
//			useSurfDes(imgs2.fc[j], imgs1.fc[i], k1, k2, cp1, cp2, g);
//			t2.stop();
//			time = t2.getTimeSec();
//			printf("GetSurfDes: %lfs\n", time);
//			count += 1;
//			if (cp1.size()>5) {
//				t2.start();
//				amh = findHomography(cp1, cp2, CV_RANSAC, 4, match_mask);
//				filterGM(g, gm, match_mask);
//				filterP(cp1, cp2, cf1, cf2, match_mask);
//				t2.stop();
//				std::printf("%d/%d FilterPointSets: %lfs\n", cf1.size(), cp1.size(), t2.getTimeSec() - time);
//
//				time = t2.getTimeSec();
//				t2.start();
//				am = estimateRigidTransform(cf1, cf2, false);
//				t2.stop();
//				printf("EstimateRigidTransform: %lfs\n", t2.getTimeSec() - time);
//				time = t2.getTimeSec();
//				PSE pe;
//				t2.start();
//				PointSetEva(cf1, cf2, pe);
//				t2.stop();
//				printf("PointSetsEva: %lfs\n", t2.getTimeSec() - time);
//				cout << am << endl << endl;
//				waitKey();
//#ifdef show
//				drawMatches(imgs2.fc[j], k1, imgs1.fc[i], k2, gm, ima, Scalar::all(-1), Scalar::all(-1),
//					vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
//#ifdef saveR
//				sprintf(tn, "%s\\%03d_%s_%s.jpg",
//					name, count, n1, n2);
//				imwrite(tn, ima);
//#else
//				resize(ima, ima, Size(1920, 540));
//				imshow("Good Matches", ima);
//				waitKey();
//#endif
//#endif // show
//				if (am.cols < 3 || am.rows < 2) {
//					std::printf("%d %s %s %d %d failed to find am\n\n",
//						count, n1, n2, i, j);
//				}
//				else {
//					estRXY e;
//					calEstR(e.rotation, am);
//					xcos = am.at<double>(0, 0);
//					xsin = am.at<double>(0, 1);
//					ax = am.at<double>(0, 2);
//					ay = am.at<double>(1, 2);
//					sqr = sqrt(xcos*xcos + xsin*xsin);
//					fprintf(f, "%s %s %lf %lf %lf %d %d %lf %lf %lf %lf %lf 21 7 %lf %lf %lf %lf %lf\n",
//						n1, n2,
//						e.rotation - (atoi(n2) - atoi(n1)),
//						e.x - 21, e.y - 7, cf1.size(), cp1.size(),
//						e.rotation, e.x, e.y, pe.mean, pe.variance, ax, ay, xcos / sqr, xsin / sqr, sqr);
//					std::printf("%d %s %s %d %d %lf %lf %lf %lf %lf\n",
//						count, n1, n2,
//						i, j, e.rotation, e.x, e.y, pe.mean, pe.variance);
//				}
//			}
//			else {
//				std::printf("%d %s %s %d %d failed to find am due to insufficient keypoints\n\n",
//					count, n1, n2, i, j);
//			}
//		}
//	}
//	t1.stop();
//	std::printf("\nTime: %lf\ns/p: %lf\n", t1.getTimeSec(), t1.getTimeSec() / (size1*size2));
//	std::fclose(f);
//}