/*
 * SL_BundleHelper.cpp
 *
 *  Created on: 2010-11-19
 *      Author: Danping Zou
 *		E-mail: dannis.zou@gmail.com
 */
#include "SL_BundleHelper.h"
#include "SL_error.h"
#include "math/SL_Matrix.h"
#include "math/SL_LinAlg.h"

#include "geometry/SL_Distortion.h"
#include "geometry/SL_Quaternion.h"
#include "extern/sba-1.6/sba.h"

#include "tools/SL_WriteRead.h"
#include "tools/SL_Debug.h"

#include <cassert>

/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
static void calcImgProjFullR(double a[5] , double qr0[4] , double t[3] , double M[3] , double n[2]) {
	double t1;
	double t11;
	double t13;
	double t17;
	double t2;
	double t22;
	double t27;
	double t3;
	double t38;
	double t46;
	double t49;
	double t5;
	double t6;
	double t8;
	double t9;
	{
		t1 = a[0];
		t2 = qr0[1];
		t3 = M[0];
		t5 = qr0[2];
		t6 = M[1];
		t8 = qr0[3];
		t9 = M[2];
		t11 = -t3 * t2 - t5 * t6 - t8 * t9;
		t13 = qr0[0];
		t17 = t13 * t3 + t5 * t9 - t8 * t6;
		t22 = t6 * t13 + t8 * t3 - t9 * t2;
		t27 = t13 * t9 + t6 * t2 - t5 * t3;
		t38 = -t5 * t11 + t13 * t22 - t27 * t2 + t8 * t17 + t[1];
		t46 = -t11 * t8 + t13 * t27 - t5 * t17 + t2 * t22 + t[2];
		t49 = 1 / t46;
		n[0] = (t1 * (-t2 * t11 + t13 * t17 - t22 * t8 + t5 * t27 + t[0]) + a[4] * t38 + a[1] * t46) * t49;
		n[1] = (t1 * a[3] * t38 + a[2] * t46) * t49;
		return;
	}
}
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
static void calcImgProjJacRTS(
		double a[5] ,
		double qr0[4] ,
		double v[3] ,
		double t[3] ,
		double M[3] ,
		double jacmRT[2][6] ,
		double jacmS[2][3]) {
	double t1;
	double t10;
	double t107;
	double t109;
	double t11;
	double t118;
	double t12;
	double t126;
	double t127;
	double t14;
	double t141;
	double t145;
	double t146;
	double t147;
	double t15;
	double t150;
	double t152;
	double t159;
	double t16;
	double t162;
	double t165;
	double t168;
	double t170;
	double t172;
	double t175;
	double t18;
	double t180;
	double t185;
	double t187;
	double t19;
	double t192;
	double t194;
	double t2;
	double t206;
	double t21;
	double t216;
	double t22;
	double t227;
	double t23;
	double t230;
	double t233;
	double t235;
	double t237;
	double t240;
	double t245;
	double t25;
	double t250;
	double t252;
	double t257;
	double t259;
	double t27;
	double t271;
	double t28;
	double t281;
	double t293;
	double t294;
	double t296;
	double t299;
	double t3;
	double t30;
	double t302;
	double t303;
	double t305;
	double t306;
	double t309;
	double t324;
	double t325;
	double t327;
	double t330;
	double t331;
	double t347;
	double t35;
	double t350;
	double t37;
	double t4;
	double t43;
	double t49;
	double t5;
	double t51;
	double t52;
	double t54;
	double t56;
	double t6;
	double t61;
	double t65;
	double t7;
	double t70;
	double t75;
	double t76;
	double t81;
	double t82;
	double t87;
	double t88;
	double t9;
	double t93;
	double t94;
	double t98;
	{
		t1 = a[0];
		t2 = v[0];
		t3 = t2 * t2;
		t4 = v[1];
		t5 = t4 * t4;
		t6 = v[2];
		t7 = t6 * t6;
		t9 = sqrt(1.0 - t3 - t5 - t7);
		t10 = 1 / t9;
		t11 = qr0[1];
		t12 = t11 * t10;
		t14 = qr0[0];
		t15 = -t12 * t2 + t14;
		t16 = M[0];
		t18 = qr0[2];
		t19 = t18 * t10;
		t21 = qr0[3];
		t22 = -t19 * t2 - t21;
		t23 = M[1];
		t25 = t10 * t21;
		t27 = -t25 * t2 + t18;
		t28 = M[2];
		t30 = -t15 * t16 - t22 * t23 - t27 * t28;
		t35 = -t9 * t11 - t2 * t14 - t4 * t21 + t6 * t18;
		t37 = -t35;
		t43 = t9 * t18 + t4 * t14 + t6 * t11 - t2 * t21;
		t49 = t9 * t21 + t6 * t14 + t2 * t18 - t11 * t4;
		t51 = -t37 * t16 - t43 * t23 - t49 * t28;
		t52 = -t15;
		t54 = t10 * t14;
		t56 = -t54 * t2 - t11;
		t61 = t9 * t14 - t2 * t11 - t4 * t18 - t6 * t21;
		t65 = t61 * t16 + t43 * t28 - t23 * t49;
		t70 = t56 * t16 + t22 * t28 - t23 * t27;
		t75 = t56 * t23 + t27 * t16 - t28 * t15;
		t76 = -t49;
		t81 = t61 * t23 + t49 * t16 - t37 * t28;
		t82 = -t27;
		t87 = t56 * t28 + t23 * t15 - t22 * t16;
		t88 = -t43;
		t93 = t61 * t28 + t37 * t23 - t43 * t16;
		t94 = -t22;
		t98 = a[4];
		t107 = t30 * t88 + t94 * t51 + t56 * t81 + t61 * t75 + t87 * t35 + t93 * t52 - t70 * t76 - t82 * t65;
		t109 = a[1];
		t118 = t30 * t76 + t82 * t51 + t56 * t93 + t61 * t87 + t70 * t88 + t65 * t94 - t35 * t75 - t81 * t52;
		t126 = t76 * t51 + t61 * t93 + t65 * t88 - t81 * t35 + t[2];
		t127 = 1 / t126;
		t141 = t51 * t88 + t61 * t81 + t93 * t35 - t65 * t76 + t[1];
		t145 = t126 * t126;
		t146 = 1 / t145;
		t147 = (t1 * (t35 * t51 + t61 * t65 + t81 * t76 - t93 * t88 + t[0]) + t98 * t141 + t126 * t109) * t146;
		jacmRT[0][0] = (t1 * (t30 * t35 + t52 * t51 + t56 * t65 + t61 * t70 + t76 * t75 + t81 * t82 - t88 * t87 - t93
				* t94) + t98 * t107 + t109 * t118) * t127 - t118 * t147;
		t150 = t1 * a[3];
		t152 = a[2];
		t159 = (t150 * t141 + t126 * t152) * t146;
		jacmRT[1][0] = (t107 * t150 + t152 * t118) * t127 - t159 * t118;
		t162 = -t12 * t4 + t21;
		t165 = -t19 * t4 + t14;
		t168 = -t25 * t4 - t11;
		t170 = -t162 * t16 - t165 * t23 - t168 * t28;
		t172 = -t162;
		t175 = -t54 * t4 - t18;
		t180 = t175 * t16 + t165 * t28 - t168 * t23;
		t185 = t175 * t23 + t168 * t16 - t162 * t28;
		t187 = -t168;
		t192 = t175 * t28 + t162 * t23 - t165 * t16;
		t194 = -t165;
		t206 = t170 * t88 + t51 * t194 + t175 * t81 + t61 * t185 + t192 * t35 + t93 * t172 - t76 * t180 - t65 * t187;
		t216 = t170 * t76 + t51 * t187 + t93 * t175 + t61 * t192 + t180 * t88 + t65 * t194 - t185 * t35 - t81 * t172;
		jacmRT[0][1] = (t1 * (t170 * t35 + t172 * t51 + t175 * t65 + t180 * t61 + t185 * t76 + t81 * t187 - t192 * t88
				- t93 * t194) + t98 * t206 + t109 * t216) * t127 - t147 * t216;
		jacmRT[1][1] = (t150 * t206 + t152 * t216) * t127 - t159 * t216;
		t227 = -t12 * t6 - t18;
		t230 = -t19 * t6 + t11;
		t233 = -t25 * t6 + t14;
		t235 = -t227 * t16 - t23 * t230 - t233 * t28;
		t237 = -t227;
		t240 = -t54 * t6 - t21;
		t245 = t240 * t16 + t230 * t28 - t233 * t23;
		t250 = t23 * t240 + t233 * t16 - t227 * t28;
		t252 = -t233;
		t257 = t240 * t28 + t227 * t23 - t230 * t16;
		t259 = -t230;
		t271 = t235 * t88 + t51 * t259 + t81 * t240 + t61 * t250 + t257 * t35 + t93 * t237 - t245 * t76 - t65 * t252;
		t281 = t235 * t76 + t51 * t252 + t240 * t93 + t61 * t257 + t245 * t88 + t259 * t65 - t250 * t35 - t81 * t237;
		jacmRT[0][2] = (t1 * (t235 * t35 + t237 * t51 + t240 * t65 + t61 * t245 + t250 * t76 + t81 * t252 - t257 * t88
				- t93 * t259) + t271 * t98 + t281 * t109) * t127 - t147 * t281;
		jacmRT[1][2] = (t150 * t271 + t281 * t152) * t127 - t159 * t281;
		jacmRT[0][3] = t127 * t1;
		jacmRT[1][3] = 0.0;
		jacmRT[0][4] = t98 * t127;
		jacmRT[1][4] = t150 * t127;
		jacmRT[0][5] = t109 * t127 - t147;
		jacmRT[1][5] = t152 * t127 - t159;
		t293 = t35 * t35;
		t294 = t61 * t61;
		t296 = t88 * t88;
		t299 = t35 * t88;
		t302 = t61 * t76;
		t303 = 2.0 * t299 + t61 * t49 - t302;
		t305 = t35 * t76;
		t306 = t61 * t88;
		t309 = t305 + 2.0 * t306 - t49 * t35;
		jacmS[0][0] = (t1 * (t293 + t294 + t49 * t76 - t296) + t98 * t303 + t109 * t309) * t127 - t147 * t309;
		jacmS[1][0] = (t150 * t303 + t152 * t309) * t127 - t159 * t309;
		t324 = t76 * t76;
		t325 = t296 + t294 + t35 * t37 - t324;
		t327 = t76 * t88;
		t330 = t61 * t35;
		t331 = 2.0 * t327 + t61 * t37 - t330;
		jacmS[0][1] = (t1 * (t299 + 2.0 * t302 - t37 * t88) + t98 * t325 + t109 * t331) * t127 - t147 * t331;
		jacmS[1][1] = (t150 * t325 + t152 * t331) * t127 - t159 * t331;
		t347 = t327 + 2.0 * t330 - t43 * t76;
		t350 = t324 + t294 + t43 * t88 - t293;
		jacmS[0][2] = (t1 * (2.0 * t305 + t61 * t43 - t306) + t98 * t347 + t350 * t109) * t127 - t147 * t350;
		jacmS[1][2] = (t150 * t347 + t152 * t350) * t127 - t159 * t350;
		return;
	}
}
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
static void calcImgProjJacRT(
		double a[5] ,
		double qr0[4] ,
		double v[3] ,
		double t[3] ,
		double M[3] ,
		double jacmRT[2][6]) {
	double t1;
	double t10;
	double t107;
	double t109;
	double t11;
	double t118;
	double t12;
	double t126;
	double t127;
	double t14;
	double t141;
	double t145;
	double t146;
	double t147;
	double t15;
	double t150;
	double t152;
	double t159;
	double t16;
	double t162;
	double t165;
	double t168;
	double t170;
	double t172;
	double t175;
	double t18;
	double t180;
	double t185;
	double t187;
	double t19;
	double t192;
	double t194;
	double t2;
	double t206;
	double t21;
	double t216;
	double t22;
	double t227;
	double t23;
	double t230;
	double t233;
	double t235;
	double t237;
	double t240;
	double t245;
	double t25;
	double t250;
	double t252;
	double t257;
	double t259;
	double t27;
	double t271;
	double t28;
	double t281;
	double t3;
	double t30;
	double t35;
	double t37;
	double t4;
	double t43;
	double t49;
	double t5;
	double t51;
	double t52;
	double t54;
	double t56;
	double t6;
	double t61;
	double t65;
	double t7;
	double t70;
	double t75;
	double t76;
	double t81;
	double t82;
	double t87;
	double t88;
	double t9;
	double t93;
	double t94;
	double t98;
	{
		t1 = a[0];
		t2 = v[0];
		t3 = t2 * t2;
		t4 = v[1];
		t5 = t4 * t4;
		t6 = v[2];
		t7 = t6 * t6;
		t9 = sqrt(1.0 - t3 - t5 - t7);
		t10 = 1 / t9;
		t11 = qr0[1];
		t12 = t11 * t10;
		t14 = qr0[0];
		t15 = -t12 * t2 + t14;
		t16 = M[0];
		t18 = qr0[2];
		t19 = t18 * t10;
		t21 = qr0[3];
		t22 = -t19 * t2 - t21;
		t23 = M[1];
		t25 = t10 * t21;
		t27 = -t25 * t2 + t18;
		t28 = M[2];
		t30 = -t15 * t16 - t22 * t23 - t27 * t28;
		t35 = -t9 * t11 - t2 * t14 - t4 * t21 + t6 * t18;
		t37 = -t35;
		t43 = t9 * t18 + t4 * t14 + t6 * t11 - t2 * t21;
		t49 = t9 * t21 + t6 * t14 + t2 * t18 - t11 * t4;
		t51 = -t37 * t16 - t43 * t23 - t49 * t28;
		t52 = -t15;
		t54 = t10 * t14;
		t56 = -t54 * t2 - t11;
		t61 = t9 * t14 - t2 * t11 - t4 * t18 - t6 * t21;
		t65 = t61 * t16 + t43 * t28 - t23 * t49;
		t70 = t56 * t16 + t22 * t28 - t23 * t27;
		t75 = t56 * t23 + t27 * t16 - t28 * t15;
		t76 = -t49;
		t81 = t61 * t23 + t49 * t16 - t37 * t28;
		t82 = -t27;
		t87 = t56 * t28 + t23 * t15 - t22 * t16;
		t88 = -t43;
		t93 = t61 * t28 + t37 * t23 - t43 * t16;
		t94 = -t22;
		t98 = a[4];
		t107 = t30 * t88 + t94 * t51 + t56 * t81 + t61 * t75 + t87 * t35 + t93 * t52 - t70 * t76 - t82 * t65;
		t109 = a[1];
		t118 = t30 * t76 + t82 * t51 + t56 * t93 + t61 * t87 + t70 * t88 + t65 * t94 - t35 * t75 - t81 * t52;
		t126 = t76 * t51 + t61 * t93 + t65 * t88 - t81 * t35 + t[2];
		t127 = 1 / t126;
		t141 = t51 * t88 + t61 * t81 + t93 * t35 - t65 * t76 + t[1];
		t145 = t126 * t126;
		t146 = 1 / t145;
		t147 = (t1 * (t35 * t51 + t61 * t65 + t81 * t76 - t93 * t88 + t[0]) + t98 * t141 + t126 * t109) * t146;
		jacmRT[0][0] = (t1 * (t30 * t35 + t52 * t51 + t56 * t65 + t61 * t70 + t76 * t75 + t81 * t82 - t88 * t87 - t93
				* t94) + t98 * t107 + t109 * t118) * t127 - t118 * t147;
		t150 = t1 * a[3];
		t152 = a[2];
		t159 = (t150 * t141 + t126 * t152) * t146;
		jacmRT[1][0] = (t107 * t150 + t152 * t118) * t127 - t159 * t118;
		t162 = -t12 * t4 + t21;
		t165 = -t19 * t4 + t14;
		t168 = -t25 * t4 - t11;
		t170 = -t162 * t16 - t165 * t23 - t168 * t28;
		t172 = -t162;
		t175 = -t54 * t4 - t18;
		t180 = t175 * t16 + t165 * t28 - t168 * t23;
		t185 = t175 * t23 + t168 * t16 - t162 * t28;
		t187 = -t168;
		t192 = t175 * t28 + t162 * t23 - t165 * t16;
		t194 = -t165;
		t206 = t170 * t88 + t51 * t194 + t175 * t81 + t61 * t185 + t192 * t35 + t93 * t172 - t76 * t180 - t65 * t187;
		t216 = t170 * t76 + t51 * t187 + t93 * t175 + t61 * t192 + t180 * t88 + t65 * t194 - t185 * t35 - t81 * t172;
		jacmRT[0][1] = (t1 * (t170 * t35 + t172 * t51 + t175 * t65 + t180 * t61 + t185 * t76 + t81 * t187 - t192 * t88
				- t93 * t194) + t98 * t206 + t109 * t216) * t127 - t147 * t216;
		jacmRT[1][1] = (t150 * t206 + t152 * t216) * t127 - t159 * t216;
		t227 = -t12 * t6 - t18;
		t230 = -t19 * t6 + t11;
		t233 = -t25 * t6 + t14;
		t235 = -t227 * t16 - t23 * t230 - t233 * t28;
		t237 = -t227;
		t240 = -t54 * t6 - t21;
		t245 = t240 * t16 + t230 * t28 - t233 * t23;
		t250 = t23 * t240 + t233 * t16 - t227 * t28;
		t252 = -t233;
		t257 = t240 * t28 + t227 * t23 - t230 * t16;
		t259 = -t230;
		t271 = t235 * t88 + t51 * t259 + t81 * t240 + t61 * t250 + t257 * t35 + t93 * t237 - t245 * t76 - t65 * t252;
		t281 = t235 * t76 + t51 * t252 + t240 * t93 + t61 * t257 + t245 * t88 + t259 * t65 - t250 * t35 - t81 * t237;
		jacmRT[0][2] = (t1 * (t235 * t35 + t237 * t51 + t240 * t65 + t61 * t245 + t250 * t76 + t81 * t252 - t257 * t88
				- t93 * t259) + t271 * t98 + t281 * t109) * t127 - t147 * t281;
		jacmRT[1][2] = (t150 * t271 + t281 * t152) * t127 - t159 * t281;
		jacmRT[0][3] = t127 * t1;
		jacmRT[1][3] = 0.0;
		jacmRT[0][4] = t98 * t127;
		jacmRT[1][4] = t150 * t127;
		jacmRT[0][5] = t109 * t127 - t147;
		jacmRT[1][5] = t152 * t127 - t159;
		return;
	}
}
#include <math.h>
static void calcImgProjJacKRTS(
		double a[5] ,
		double qr0[4] ,
		double v[3] ,
		double t[3] ,
		double M[3] ,
		double jacmKRT[2][11] ,
		double jacmS[2][3]) {
	double t1;
	double t102;
	double t107;
	double t109;
	double t11;
	double t114;
	double t116;
	double t120;
	double t129;
	double t13;
	double t131;
	double t140;
	double t148;
	double t149;
	double t15;
	double t150;
	double t152;
	double t154;
	double t161;
	double t164;
	double t167;
	double t17;
	double t170;
	double t172;
	double t174;
	double t177;
	double t18;
	double t182;
	double t187;
	double t189;
	double t194;
	double t196;
	double t2;
	double t208;
	double t218;
	double t229;
	double t232;
	double t235;
	double t237;
	double t239;
	double t24;
	double t242;
	double t247;
	double t25;
	double t252;
	double t254;
	double t259;
	double t261;
	double t273;
	double t283;
	double t295;
	double t296;
	double t298;
	double t3;
	double t301;
	double t304;
	double t305;
	double t307;
	double t308;
	double t31;
	double t311;
	double t32;
	double t326;
	double t327;
	double t329;
	double t332;
	double t333;
	double t34;
	double t349;
	double t35;
	double t352;
	double t4;
	double t41;
	double t45;
	double t5;
	double t50;
	double t51;
	double t56;
	double t57;
	double t6;
	double t60;
	double t66;
	double t67;
	double t68;
	double t74;
	double t76;
	double t78;
	double t79;
	double t8;
	double t81;
	double t83;
	double t85;
	double t87;
	double t89;
	double t9;
	double t91;
	double t93;
	double t95;
	double t97;
	{
		t1 = v[0];
		t2 = t1 * t1;
		t3 = v[1];
		t4 = t3 * t3;
		t5 = v[2];
		t6 = t5 * t5;
		t8 = sqrt(1.0 - t2 - t4 - t6);
		t9 = qr0[1];
		t11 = qr0[0];
		t13 = qr0[3];
		t15 = qr0[2];
		t17 = t8 * t9 + t11 * t1 + t13 * t3 - t5 * t15;
		t18 = M[0];
		t24 = t8 * t15 + t3 * t11 + t5 * t9 - t13 * t1;
		t25 = M[1];
		t31 = t8 * t13 + t5 * t11 + t1 * t15 - t3 * t9;
		t32 = M[2];
		t34 = -t17 * t18 - t24 * t25 - t31 * t32;
		t35 = -t17;
		t41 = t11 * t8 - t1 * t9 - t3 * t15 - t5 * t13;
		t45 = t41 * t18 + t24 * t32 - t31 * t25;
		t50 = t41 * t25 + t31 * t18 - t17 * t32;
		t51 = -t31;
		t56 = t41 * t32 + t17 * t25 - t24 * t18;
		t57 = -t24;
		t60 = t34 * t35 + t41 * t45 + t50 * t51 - t56 * t57 + t[0];
		t66 = t34 * t51 + t41 * t56 + t45 * t57 - t50 * t35 + t[2];
		t67 = 1 / t66;
		jacmKRT[0][0] = t60 * t67;
		t68 = a[3];
		t74 = t34 * t57 + t41 * t50 + t56 * t35 - t45 * t51 + t[1];
		jacmKRT[1][0] = t68 * t74 * t67;
		jacmKRT[0][1] = 1.0;
		jacmKRT[1][1] = 0.0;
		jacmKRT[0][2] = 0.0;
		jacmKRT[1][2] = 1.0;
		jacmKRT[0][3] = 0.0;
		t76 = a[0];
		jacmKRT[1][3] = t76 * t74 * t67;
		jacmKRT[0][4] = t74 * t67;
		jacmKRT[1][4] = 0.0;
		t78 = 1 / t8;
		t79 = t78 * t9;
		t81 = -t79 * t1 + t11;
		t83 = t78 * t15;
		t85 = -t83 * t1 - t13;
		t87 = t78 * t13;
		t89 = -t87 * t1 + t15;
		t91 = -t81 * t18 - t85 * t25 - t89 * t32;
		t93 = -t81;
		t95 = t78 * t11;
		t97 = -t95 * t1 - t9;
		t102 = t97 * t18 + t85 * t32 - t89 * t25;
		t107 = t97 * t25 + t89 * t18 - t81 * t32;
		t109 = -t89;
		t114 = t97 * t32 + t81 * t25 - t85 * t18;
		t116 = -t85;
		t120 = a[4];
		t129 = t91 * t57 + t34 * t116 + t97 * t50 + t41 * t107 + t114 * t35 + t56 * t93 - t102 * t51 - t45 * t109;
		t131 = a[1];
		t140 = t91 * t51 + t34 * t109 + t97 * t56 + t41 * t114 + t102 * t57 + t45 * t116 - t107 * t35 - t50 * t93;
		t148 = t66 * t66;
		t149 = 1 / t148;
		t150 = (t76 * t60 + t120 * t74 + t131 * t66) * t149;
		jacmKRT[0][5] = (t76 * (t91 * t35 + t34 * t93 + t97 * t45 + t41 * t102 + t107 * t51 + t50 * t109 - t114 * t57
				- t56 * t116) + t129 * t120 + t131 * t140) * t67 - t150 * t140;
		t152 = t76 * t68;
		t154 = a[2];
		t161 = (t152 * t74 + t154 * t66) * t149;
		jacmKRT[1][5] = (t152 * t129 + t154 * t140) * t67 - t161 * t140;
		t164 = -t79 * t3 + t13;
		t167 = -t83 * t3 + t11;
		t170 = -t87 * t3 - t9;
		t172 = -t164 * t18 - t167 * t25 - t170 * t32;
		t174 = -t164;
		t177 = -t95 * t3 - t15;
		t182 = t177 * t18 + t167 * t32 - t170 * t25;
		t187 = t177 * t25 + t170 * t18 - t164 * t32;
		t189 = -t170;
		t194 = t177 * t32 + t164 * t25 - t167 * t18;
		t196 = -t167;
		t208 = t172 * t57 + t34 * t196 + t177 * t50 + t41 * t187 + t194 * t35 + t56 * t174 - t182 * t51 - t45 * t189;
		t218 = t172 * t51 + t34 * t189 + t177 * t56 + t41 * t194 + t182 * t57 + t45 * t196 - t187 * t35 - t50 * t174;
		jacmKRT[0][6] = (t76 * (t172 * t35 + t34 * t174 + t177 * t45 + t41 * t182 + t187 * t51 + t50 * t189 - t194
				* t57 - t56 * t196) + t120 * t208 + t131 * t218) * t67 - t150 * t218;
		jacmKRT[1][6] = (t152 * t208 + t154 * t218) * t67 - t161 * t218;
		t229 = -t79 * t5 - t15;
		t232 = -t83 * t5 + t9;
		t235 = -t87 * t5 + t11;
		t237 = -t229 * t18 - t232 * t25 - t235 * t32;
		t239 = -t229;
		t242 = -t95 * t5 - t13;
		t247 = t242 * t18 + t232 * t32 - t235 * t25;
		t252 = t242 * t25 + t235 * t18 - t229 * t32;
		t254 = -t235;
		t259 = t242 * t32 + t229 * t25 - t232 * t18;
		t261 = -t232;
		t273 = t237 * t57 + t261 * t34 + t242 * t50 + t41 * t252 + t259 * t35 + t56 * t239 - t247 * t51 - t45 * t254;
		t283 = t237 * t51 + t34 * t254 + t242 * t56 + t41 * t259 + t247 * t57 + t45 * t261 - t252 * t35 - t50 * t239;
		jacmKRT[0][7] = (t76 * (t237 * t35 + t34 * t239 + t242 * t45 + t41 * t247 + t252 * t51 + t50 * t254 - t259
				* t57 - t56 * t261) + t120 * t273 + t131 * t283) * t67 - t150 * t283;
		jacmKRT[1][7] = (t152 * t273 + t154 * t283) * t67 - t161 * t283;
		jacmKRT[0][8] = t76 * t67;
		jacmKRT[1][8] = 0.0;
		jacmKRT[0][9] = t120 * t67;
		jacmKRT[1][9] = t152 * t67;
		jacmKRT[0][10] = t131 * t67 - t150;
		jacmKRT[1][10] = t154 * t67 - t161;
		t295 = t35 * t35;
		t296 = t41 * t41;
		t298 = t57 * t57;
		t301 = t35 * t57;
		t304 = t41 * t51;
		t305 = 2.0 * t301 + t41 * t31 - t304;
		t307 = t35 * t51;
		t308 = t41 * t57;
		t311 = t307 + 2.0 * t308 - t31 * t35;
		jacmS[0][0] = (t76 * (t295 + t296 + t31 * t51 - t298) + t120 * t305 + t131 * t311) * t67 - t150 * t311;
		jacmS[1][0] = (t152 * t305 + t154 * t311) * t67 - t161 * t311;
		t326 = t51 * t51;
		t327 = t298 + t296 + t17 * t35 - t326;
		t329 = t57 * t51;
		t332 = t41 * t35;
		t333 = 2.0 * t329 + t41 * t17 - t332;
		jacmS[0][1] = (t76 * (t301 + 2.0 * t304 - t17 * t57) + t120 * t327 + t131 * t333) * t67 - t150 * t333;
		jacmS[1][1] = (t152 * t327 + t154 * t333) * t67 - t161 * t333;
		t349 = t329 + 2.0 * t332 - t24 * t51;
		t352 = t326 + t296 + t24 * t57 - t295;
		jacmS[0][2] = (t76 * (2.0 * t307 + t41 * t24 - t308) + t120 * t349 + t131 * t352) * t67 - t150 * t352;
		jacmS[1][2] = (t152 * t349 + t154 * t352) * t67 - t161 * t352;
		return;
	}
}

/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
void img_projsRTS_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata) {
	register int i, j;
	int cnp, pnp, mnp;
	double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *Kparms, *pr0;
	//int n;
	int m, nnz, Asz, Bsz, ABsz;
	sbaGlobs *gl;

	gl = (sbaGlobs *) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;
	Kparms = gl->intrcalib;

	m = idxij->nc;
	pa = p;
	pb = p + m * cnp;
	Asz = mnp * cnp;
	Bsz = mnp * pnp;
	ABsz = Asz + Bsz;

	for (j = 0; j < m; ++j) {
		/* idx2-th camera parameters */
		pqr = pa + j * cnp;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate
		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, idx1=0...n-1 */
		for (i = 0; i < nnz; ++i) {
			ppt = pb + rcsubs[i] * pnp;
			pA = jac + idxij->val[rcidxs[i]] * ABsz; // set pA to point to A_ij
			pB = pA + Asz; // set pB to point to B_ij
			calcImgProjJacRTS(Kparms, pr0, pqr, pt, ppt, (double(*)[6]) pA, (double(*)[3]) pB); // evaluate dQ/da, dQ/db in pA, pB
		}
	}
}
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
void img_projsRTS_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata) {
	register int i, j;
	int cnp, pnp, mnp;
	double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *Kparms, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
	//int n;
	int m, nnz;
	sbaGlobs *gl;

	gl = (sbaGlobs*) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;
	Kparms = gl->intrcalib;

	m = idxij->nc;
	pa = p;
	pb = p + m * cnp;

	for (j = 0; j < m; ++j) {
		/* idx2-th camera parameters */
		pqr = pa + j * cnp;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate
		_MK_QUAT_FRM_VEC(lrot, pqr);
		quatMultFast(lrot, pr0, trot); // trot=lrot*pr0
		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, idx1=0...n-1 */

		for (i = 0; i < nnz; ++i) {
			ppt = pb + rcsubs[i] * pnp;
			pmeas = hx + idxij->val[rcidxs[i]] * mnp; // set pmeas to point to hx_ij
			calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
		}
	}
}
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
void img_projsRT_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata) {
	register int i, j;
	int cnp, pnp, mnp;
	double *pqr, *pt, *ppt, *pmeas, *Kparms, *ptparams, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
	//int n;
	int m, nnz;
	sbaGlobs *gl;

	gl = (sbaGlobs *) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;
	Kparms = gl->intrcalib;
	ptparams = gl->ptparams;

	m = idxij->nc;

	for (j = 0; j < m; ++j) {
		/* idx2-th camera parameters */
		pqr = p + j * cnp;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate
		_MK_QUAT_FRM_VEC(lrot, pqr);
		quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, idx1=0...n-1 */

		for (i = 0; i < nnz; ++i) {
			ppt = ptparams + rcsubs[i] * pnp;
			pmeas = hx + idxij->val[rcidxs[i]] * mnp; // set pmeas to point to hx_ij
			calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
		}
	}
}
/*
 * copy from 'eucsbademo.c' in sba-1.6.1 
 */
void img_projsRT_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata) {
	register int i, j;
	int cnp, pnp, mnp;
	double *pqr, *pt, *ppt, *pA, *Kparms, *ptparams, *pr0;
	//int n;
	int m, nnz, Asz;
	sbaGlobs *gl;

	gl = (sbaGlobs *) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;
	Kparms = gl->intrcalib;
	ptparams = gl->ptparams;

	//n=idxij->nr;
	m = idxij->nc;
	Asz = mnp * cnp;

	for (j = 0; j < m; ++j) {
		/* idx2-th camera parameters */
		pqr = p + j * cnp;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate

		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, idx1=0...n-1 */

		for (i = 0; i < nnz; ++i) {
			ppt = ptparams + rcsubs[i] * pnp;
			pA = jac + idxij->val[rcidxs[i]] * Asz; // set pA to point to A_ij

			calcImgProjJacRT(Kparms, pr0, pqr, pt, ppt, (double(*)[6]) pA); // evaluate dQ/da in pA
		}
	}
}
/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_projsKRTS_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *hx , void *adata) {
	register int i, j;
	int cnp, pnp, mnp;
	double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *pcalib, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
	//int n;
	int m, nnz;
	sbaGlobs *gl;

	gl = (sbaGlobs *) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;

	//n=idxij->nr;
	m = idxij->nc;
	pa = p;
	pb = p + m * cnp;

	for (j = 0; j < m; ++j) {
		/* j-th camera parameters */
		pcalib = pa + j * cnp;
		pqr = pcalib + 5;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate
		_MK_QUAT_FRM_VEC(lrot, pqr);
		quatMultFast(lrot, pr0, trot); // trot=lrot*pr0

		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

		for (i = 0; i < nnz; ++i) {
			ppt = pb + rcsubs[i] * pnp;
			pmeas = hx + idxij->val[rcidxs[i]] * mnp; // set pmeas to point to hx_ij

			calcImgProjFullR(pcalib, trot, pt, ppt, pmeas); // evaluate Q in pmeas
			//calcImgProj(pcalib, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
		}
	}
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
void img_projsKRTS_jac_x(double *p , struct sba_crsm *idxij , int *rcidxs , int *rcsubs , double *jac , void *adata) {
	register int i, j, ii, jj;
	int cnp, pnp, mnp, ncK;
	double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *pcalib, *pr0;
	//int n;
	int m, nnz, Asz, Bsz, ABsz;
	sbaGlobs *gl;

	gl = (sbaGlobs*) adata;
	cnp = gl->cnp;
	pnp = gl->pnp;
	mnp = gl->mnp;
	ncK = gl->nccalib;

	//n=idxij->nr;
	m = idxij->nc;
	pa = p;
	pb = p + m * cnp;
	Asz = mnp * cnp;
	Bsz = mnp * pnp;
	ABsz = Asz + Bsz;

	for (j = 0; j < m; ++j) {
		/* j-th camera parameters */
		pcalib = pa + j * cnp;
		pqr = pcalib + 5;
		pt = pqr + 3; // quaternion vector part has 3 elements
		pr0 = gl->rot0params + j * FULLQUATSZ; // full quat for initial rotation estimate

		nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

		for (i = 0; i < nnz; ++i) {
			ppt = pb + rcsubs[i] * pnp;
			pA = jac + idxij->val[rcidxs[i]] * ABsz; // set pA to point to A_ij
			pB = pA + Asz; // set pB to point to B_ij

			calcImgProjJacKRTS(pcalib, pr0, pqr, pt, ppt, (double(*)[5 + 6]) pA, (double(*)[3]) pB); // evaluate dQ/da, dQ/db in pA, pB

			/* clear the columns of the Jacobian corresponding to fixed calibration parameters */
			if (ncK) {
				int jj0 = 5 - ncK;

				for (ii = 0; ii < mnp; ++ii, pA += cnp)
					for (jj = jj0; jj < 5; ++jj)
						pA[jj] = 0.0; // pA[ii*cnp+jj]=0.0;
			}
		}
	}
}
void baRTS(
		int npts ,
		const double* meas ,
		const double* covx ,
		const char* vmask ,
		double* pts ,
		int nv ,
		int nvcon ,
		double* Rs ,
		double* ts ,
		int max_iter) {
	const int cnp = 6; //number of camera parameters, 3 for R, 3 for t
	const int pnp = 3; //dimension of 3D points
	const int mnp = 2; //dimension of 2D measurements

	sbaGlobs globs;
	globs.cnp = cnp;
	globs.pnp = pnp;
	globs.mnp = mnp;

	globs.rot0params = new double[FULLQUATSZ * nv];

	int i;
	for (i = 0; i < nv; ++i) {
		mat2quat(Rs + 9 * i, globs.rot0params + 4 * i);
	}

	globs.camparams = 0;
	globs.ptparams = 0;
	globs.intrcalib = new double[5];
	globs.intrcalib[0] = 1;
	globs.intrcalib[1] = 0;
	globs.intrcalib[2] = 0;
	globs.intrcalib[3] = 1;
	globs.intrcalib[4] = 0;

	double opts[SBA_OPTSSZ];

	opts[0] = SBA_INIT_MU;
	opts[1] = SBA_STOP_THRESH;
	opts[2] = SBA_STOP_THRESH;
	opts[3] = SBA_STOP_THRESH;
	//opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
	opts[4] = 0.0;

	//set the parameter vector
	double* param_vec = new double[nv * cnp + npts * pnp];
	for (i = 0; i < nv; ++i) {
		param_vec[i * cnp] = 0;
		param_vec[i * cnp + 1] = 0;
		param_vec[i * cnp + 2] = 0;
		param_vec[i * cnp + 3] = ts[i * 3];
		param_vec[i * cnp + 4] = ts[i * 3 + 1];
		param_vec[i * cnp + 5] = ts[i * 3 + 2];
	}
	double* param_pts = param_vec + nv * cnp;

	memcpy(param_pts, pts, npts * 3 * sizeof(double));

	double info_data[SBA_INFOSZ];
	sba_motstr_levmar_x(
			npts,
			0,
			nv,
			nvcon,
			const_cast<char*> (vmask),
			param_vec,
			cnp,
			pnp,
			const_cast<double*> (meas),
			covx == 0 ? 0 : const_cast<double*> (covx),
			mnp,
			img_projsRTS_x,
			img_projsRTS_jac_x,
			&globs,
			max_iter,
			0,
			opts,
			info_data);

	//test
	logInfo(
			"BA_RTS: initial error:%lf, final error:%lf #iterations:%lf stop reason:%lf\n",
			info_data[0],
			info_data[1],
			info_data[5],
			info_data[6]);

	//copy back the results	
	double dq[4], q[4];
	for (i = 0; i < nv; ++i) {
		//rotation
		_MK_QUAT_FRM_VEC(dq,param_vec+i*cnp);
		quatMultFast(dq, globs.rot0params + i * FULLQUATSZ, q);
		quat2mat(q, Rs + i * 9);

		//translation
		double* pts = ts + i * 3;
		pts[0] = param_vec[i * cnp + 3];
		pts[1] = param_vec[i * cnp + 4];
		pts[2] = param_vec[i * cnp + 5];
	}

	memcpy(pts, param_pts, npts * 3 * sizeof(double));

	delete[] param_vec;
	delete[] globs.rot0params;
	delete[] globs.intrcalib;
}

//#include "SL_Debug.h"
//static double compute_reprojection_error(
//		int v_id ,
//		int nv ,
//		const double* Rs ,
//		const double* ts ,
//		int npts ,
//		const double* pts ,
//		const double* meas ,
//		const char* vmask) {
//	const double* R = Rs + v_id * 9;
//	const double* t = ts + v_id * 3;
//
//	Mat_i indices(npts, nv);
//	indices.fill(-1);
//	int k = 0;
//	for (int i = 0; i < npts * nv; i++) {
//		if (vmask[i] == 1)
//			indices.data[i] = k++;
//	}
//
//	double err = 0;
//	double m[2];
//	for (int i = 0; i < npts; i++) {
//		if (vmask[i * nv + v_id] == 0)
//			continue;
//		project(R, t, pts + 3 * i, m);
//
//		int ind = indices[i * nv + v_id];
//		const double* pmeas = meas + 2 * ind;
//		//info("(%lf,%lf) - (%lf,%lf)\n", m[0], m[1], pmeas[0], pmeas[1]);
//		err += (m[0] - pmeas[0]) * (m[0] - pmeas[0]) + (m[1] - pmeas[1]) * (m[1] - pmeas[1]);
//	}
//	return sqrt(err / npts);
//}
void baRTS(
		const double* K ,
		const double* kc ,
		const double* iK ,
		const double* k_ud ,
		int npts ,
		int nmeas ,
		const double* meas ,
		const double* covx ,
		const char* vmask ,
		double* pts ,
		int nv ,
		int nvcon ,
		double *Rs ,
		double* ts ,
		int max_iter) {

	assert(npts > 0 && nv > 0);
	double* meas_norm = new double[nmeas * 2];
	undistorNormPoints(iK, k_ud, nmeas, meas, meas_norm);

	//	info("reprojection error before BA:\n");
	//	for (int v = 0; v < nv; ++v) {
	//		draw_reprojection_points(K, kc, v, nv, Rs, ts, npts, pts, meas_norm, vmask);
	//		double err = compute_reprojection_error(v, nv, Rs, ts, npts, pts, meas_norm, vmask);
	//		info("view %d : %lf\n", v, err);
	//	}

	baRTS(npts, meas_norm, covx, vmask, pts, nv, nvcon, Rs, ts, max_iter);

	//	info("reprojection error after BA:\n");
	//	for (int v = 0; v < nv; ++v) {
	//		draw_reprojection_points(K, kc, v, nv, Rs, ts, npts, pts, meas_norm, vmask);
	//		double err = compute_reprojection_error(v, nv, Rs, ts, npts, pts, meas_norm, vmask);
	//		info("view %d : %lf\n", v, err);
	//	}
	//
	//	cv::waitKey(-1);

	delete[] meas_norm;
}

void baRTS2(
		int npts ,
		const double* m0 ,
		const double* covx0 ,
		const double* m1 ,
		const double* covx1 ,
		const double* R0 ,
		const double* t0 ,
		const double* R1 ,
		const double* t1 ,
		const double* M ,
		double* outR0 ,
		double* outT0 ,
		double* outR1 ,
		double* outT1 ,
		double* outM ,
		int matIter) {

	assert(npts>0);
	assert( (covx0 == 0 && covx1 == 0) || (covx0 != 0 && covx1 != 0));

	//set the parameters
	const int cnp = 6;
	const int pnp = 3;
	const int mnp = 2;

	sbaGlobs globs;
	globs.cnp = cnp;
	globs.pnp = pnp;
	globs.mnp = mnp;

	globs.rot0params = new double[FULLQUATSZ * 2];

	mat2quat(R0, globs.rot0params);
	mat2quat(R1, globs.rot0params + 4);

	globs.camparams = 0;
	globs.ptparams = 0;

	globs.intrcalib = new double[5];
	globs.intrcalib[0] = 1;
	globs.intrcalib[1] = 0;
	globs.intrcalib[2] = 0;
	globs.intrcalib[3] = 1;
	globs.intrcalib[4] = 0;

	double opts[SBA_OPTSSZ];

	opts[0] = SBA_INIT_MU;
	opts[1] = SBA_STOP_THRESH;
	opts[2] = SBA_STOP_THRESH;
	opts[3] = SBA_STOP_THRESH;
	//opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
	opts[4] = 0.0;

	char* vmask = new char[npts * 2];
	memset(vmask, 1, sizeof(char) * npts * 2);

	//fill initial parameters
	double* param_vec = new double[2 * cnp + npts * pnp];
	param_vec[0] = 0;
	param_vec[1] = 0;
	param_vec[2] = 0;
	param_vec[3] = t0[0];
	param_vec[4] = t0[1];
	param_vec[5] = t0[2];
	param_vec[6] = 0;
	param_vec[7] = 0;
	param_vec[8] = 0;
	param_vec[9] = t1[0];
	param_vec[10] = t1[1];
	param_vec[11] = t1[2];

	int i;
	double* ppar = param_vec + 12;
	const double* pM = M;
	for (i = 0; i < npts; ++i) {
		*ppar = *pM;
		*(ppar + 1) = *(pM + 1);
		*(ppar + 2) = *(pM + 2);
		pM += 3;
		ppar += 3;
	}

	//fill the measurement vector
	double* meas_vec = new double[npts * 2 * mnp];
	double* pmeas = meas_vec;
	const double* pm0 = m0;
	const double* pm1 = m1;
	for (i = 0; i < npts; ++i) {
		*pmeas = *pm0;
		*(pmeas + 1) = *(pm0 + 1);
		*(pmeas + 2) = *pm1;
		*(pmeas + 3) = *(pm1 + 1);
		pmeas += 4;
		pm0 += 2;
		pm1 += 2;
	}
	double info_data[SBA_INFOSZ];

	if (covx0 != 0 && covx1 != 0) {
		double* covx_vec = new double[npts * 2 * 4];
		double* pcovx = covx_vec;
		const double* pcovx0 = covx0;
		const double* pcovx1 = covx1;
		for (i = 0; i < npts; ++i) {
			pcovx[0] = *pcovx0;
			pcovx[1] = 0;
			pcovx[2] = 0;
			pcovx[3] = *pcovx0;

			pcovx[4] = *pcovx1;
			pcovx[5] = 0;
			pcovx[6] = 0;
			pcovx[7] = *pcovx1;

			pcovx += 8;
			pcovx0++;
			pcovx1++;
		}

		sba_motstr_levmar_x(
				npts,
				1,
				2,
				1,
				vmask,
				param_vec,
				cnp,
				pnp,
				meas_vec,
				covx_vec,
				mnp,
				img_projsRTS_x,
				img_projsRTS_jac_x,
				&globs,
				matIter,
				0,
				opts,
				info_data);
		delete[] covx_vec;
	} else
		sba_motstr_levmar_x(
				npts,
				1,
				2,
				1,
				vmask,
				param_vec,
				cnp,
				pnp,
				meas_vec,
				0,
				mnp,
				img_projsRTS_x,
				img_projsRTS_jac_x,
				&globs,
				matIter,
				0,
				opts,
				info_data);

	//test
	logInfo(
			"initial error:%lf, final error:%lf #iterations:%lf stop reason:%lf\n",
			info_data[0],
			info_data[1],
			info_data[5],
			info_data[6]);

	//copy back the results	
	double dq[4], q[4];
	_MK_QUAT_FRM_VEC(dq, param_vec);
	quatMultFast(dq, globs.rot0params, q);
	quat2mat(q, outR0);
	outT0[0] = param_vec[3];
	outT0[1] = param_vec[4];
	outT0[2] = param_vec[5];

	_MK_QUAT_FRM_VEC(dq, param_vec+6);
	quatMultFast(dq, globs.rot0params + 4, q);
	quat2mat(q, outR1);
	outT1[0] = param_vec[9];
	outT1[1] = param_vec[10];
	outT1[2] = param_vec[11];

	ppar = param_vec + 12;
	double *pM_ = outM;
	for (i = 0; i < npts; ++i) {
		*pM_ = *ppar;
		*(pM_ + 1) = *(ppar + 1);
		*(pM_ + 2) = *(ppar + 2);
		pM_ += 3;
		ppar += 3;
	}

	delete[] meas_vec;
	delete[] param_vec;
	delete[] vmask;
	delete[] globs.rot0params;
	delete[] globs.intrcalib;
}
void baRTS2(
		const double* iK ,
		const double* k_ud ,
		int npts ,
		const double* m0 ,
		const double* covx0 ,
		const double* m1 ,
		const double* covx1 ,
		const double* R0 ,
		const double* t0 ,
		const double* R1 ,
		const double* t1 ,
		const double* M ,
		double* R0_ ,
		double* t0_ ,
		double* R1_ ,
		double* t1_ ,
		double* M_ ,
		int max_iter) {

	double* norm_m0 = new double[npts * 2];
	double* norm_m1 = new double[npts * 2];

	undistorNormPoints(iK, k_ud, npts, m0, norm_m0);
	undistorNormPoints(iK, k_ud, npts, m1, norm_m1);

	baRTS2(npts, norm_m0, covx0, norm_m1, covx1, R0, t0, R1, t1, M, R0_, t0_, R1_, t1_, M_, max_iter);

	delete[] norm_m0;
	delete[] norm_m1;
}
