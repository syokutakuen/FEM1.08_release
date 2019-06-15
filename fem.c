/*
 * ���܂�ȒP�ł͂Ȃ��񎟌��v�f��FEM�v���O�����B
 *  ��ҁF����_�s(�Y�@�G���W�j�A�����O �}�l�W���[)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>

/*
 *  �ėp�}�N��
 */
#define ON  (1)
#define OFF (0)
#ifdef __GCC__
#define offsetof(st, m) ((size_t)&(((st *)0)->m))
#endif
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

#ifdef __LLU__
#define LLU0	"%lu"
#define LLU(x)	"%" #x "lu"
#else
#define LLU0	"%I64u"
#define LLU(x)	"%"#x"I64u"
#endif

//MODULE INTEGRAL
/*******************************************************************************
 *
 * �K�E�X�̐��l�ϕ��_�Əd�݂̒�`
 *
 ******************************************************************************/
/*
 * ���v�f�̐ϕ��_���W
 *  �e�`�󂲂Ǝ������ƂɌʂ̔z��Ƃ��Ă��邪�A
 *  ���p����Ƃ��͎��������܂��Ă���̂ŁA
 *  ���̕����������������낤�Ɣ��f�������Ƃɂ��B
 */
/* ����1�`3���̐ϕ��_�ʒu�̒�`�B */
#define __LL1_it  0.577350269189626
#define __LL2_it  0.774596669241483
#define __LL3_it1 0.861136311594053
#define __LL3_it2 0.339981043584856

/* �ϕ��_�̈ʒu�̔z��B�e�����p�Ɍʂ̔z��Ƃ��Ē�`����B */
// static double LL0_it [] = { 0, };
static double LL1_it [] = { -__LL1_it, __LL1_it, };
static double LL2_it [] = { -__LL2_it, 0, __LL2_it, };
static double LL3_it [] = { -__LL3_it1, -__LL3_it2, __LL3_it2, __LL3_it1, };

/* 3���v�f�̂��߂̏d�݂̒l�B�������Ȃ̂Ń}�N����`����B */
#define __LL3_wt1 0.347854845137454
#define __LL3_wt2 0.652145154862456

/* 0�`3���̏d�݂̔z��B�ϕ��_�ʒu���l�A�ʂ̔z��Ƃ���B */
// static double LL0_wt [] = { 2 };
static double LL1_wt [] = { 1, 1, };
static double LL2_wt [] = { 5./9, 8./9, 5./9, };
static double LL3_wt [] = { __LL3_wt1, __LL3_wt2, __LL3_wt2, __LL3_wt1, };

/* �O�p�`�v�f�̐ϕ��_�Ǝ����B�O�p�`�v�f�ɂ�0���̂��̂͂Ȃ��̂ŏȗ�����B */
static double Tr1_it [] = { 1./3., 1./3., };
static double Tr2_it [] = { 0.5, 0.5, 0, 0.5, 0.5, 0, };
static double Tr3_it [] = { 1./3., 1./3., 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 1, 0, 0, 1, 0, 0, };
static double Tr1_wt [] = { .5, };
static double Tr2_wt [] = { 1./6., 1./6., 1./6., };
static double Tr3_wt [] = { 9./40, 1./15, 1./15, 1./15, 1./40., 1./40., 1./40., };

/* �l�p�`�v�f�̐ϕ��_�ʒu�Əd�݁B���v�f�̃f�[�^�[����쐬�����B */
static double Q0_it [] = { 0, 0, };
static double Q1_it [] = {
	-__LL1_it, -__LL1_it, __LL1_it, -__LL1_it,
	-__LL1_it,  __LL1_it, __LL1_it,  __LL1_it,
};
static double Q2_it [] = {
	-__LL2_it, -__LL2_it, 0, -__LL2_it, __LL2_it, -__LL2_it,
	-__LL2_it, 0,         0, 0,         __LL2_it, 0,
	-__LL2_it,  __LL2_it, 0,  __LL2_it, __LL2_it,  __LL2_it,
};
static double Q3_it [] = {
	-__LL3_it1, -__LL3_it1, -__LL3_it2, -__LL3_it1, __LL3_it2, -__LL3_it1,  __LL3_it1, -__LL3_it1,
	-__LL3_it1, -__LL3_it2, -__LL3_it2, -__LL3_it2, __LL3_it2, -__LL3_it2,  __LL3_it1, -__LL3_it2,
	-__LL3_it1,  __LL3_it2, -__LL3_it2,  __LL3_it2, __LL3_it2,  __LL3_it2,  __LL3_it1,  __LL3_it2,
	-__LL3_it1,  __LL3_it1, -__LL3_it2,  __LL3_it1, __LL3_it2,  __LL3_it1,  __LL3_it1,  __LL3_it1,
};

static double Q0_wt [] = { 4, };
static double Q1_wt [] = { 1, 1, 1, 1, };
static double Q2_wt [] = { 25./81, 40./81, 25./81, 40./81, 64./81, 40./81, 25./81, 40./81, 25./81, };

#define __LL3_wt11 (__LL3_wt1*__LL3_wt1)
#define __LL3_wt12 (__LL3_wt1*__LL3_wt2)
#define __LL3_wt21 __LL3_wt12
#define __LL3_wt22 (__LL3_wt2*__LL3_wt2)
static double Q3_wt [] = {
	__LL3_wt11, __LL3_wt12, __LL3_wt12, __LL3_wt11,
	__LL3_wt21, __LL3_wt22, __LL3_wt22, __LL3_wt21,
	__LL3_wt21, __LL3_wt22, __LL3_wt22, __LL3_wt21,
	__LL3_wt11, __LL3_wt12, __LL3_wt12, __LL3_wt11,
};

// �Ȃ��׏d��������Ƃ��̕ӂɂ�����ϕ��_�ʒu
#define SQRT_2 1.4142135623730950488016887242097
#define CVT_PT1(x) (SQRT_2*((x)+1)/2)
#define CVT_PT2(x) (((x)+1)/2)
static double TrE1_it [] = {
//      t1                 t2                    t3
	CVT_PT1(-__LL1_it), CVT_PT1( __LL1_it), // 0
	CVT_PT1( __LL1_it), CVT_PT1(-__LL1_it), // 0
	0.0,                CVT_PT2( __LL1_it), // CVT_PT2(-__LL1_it)
	0.0,                CVT_PT2(-__LL1_it), // CVT_PT2( __LL1_it)
	CVT_PT2(-__LL1_it), 0.0,                // CVT_PT2( __LL1_it)
	CVT_PT2( __LL1_it), 0.0,                // CVT_PT2(-__LL1_it)
};
static double TrE2_it [] = {
//      t1                 t2                    t3
	CVT_PT1(-__LL2_it), CVT_PT1( __LL2_it), // 0
	0.5,                0.5,                // 0
	CVT_PT1( __LL2_it), CVT_PT1(-__LL2_it), // 0
	0.0,                CVT_PT2( __LL2_it), // CVT_PT2(-__LL2_it)
	0.5,                0,                  // 0.5
	0.0,                CVT_PT2(-__LL2_it), // CVT_PT2( __LL2_it)
	CVT_PT2(-__LL2_it), 0.0,                // CVT_PT2(-__LL2_it)
	0.0,                0.5,                // 0.5
	CVT_PT2( __LL2_it), 0.0,                // CVT_PT2( __LL2_it)
};
static double TrE3_it [] = {
//      t1                 t2                    t3
	CVT_PT1(-__LL3_it1), CVT_PT1( __LL3_it1), // 0
	CVT_PT1( __LL3_it2), CVT_PT1(-__LL3_it2), // 0
	CVT_PT1(-__LL3_it2), CVT_PT1( __LL3_it2), // 0
	CVT_PT1( __LL3_it1), CVT_PT1(-__LL3_it1), // 0
	0.0,                 CVT_PT2( __LL3_it1), // CVT_PT2(-__LL3_it1)
	0.0,                 CVT_PT2( __LL3_it2), // CVT_PT2(-__LL3_it2)
	0.0,                 CVT_PT2(-__LL3_it2), // CVT_PT2( __LL3_it2)
	0.0,                 CVT_PT2(-__LL3_it1), // CVT_PT2( __LL3_it1)
	CVT_PT2(-__LL3_it1), 0.0,                 // CVT_PT2( __LL3_it1)
	CVT_PT2(-__LL3_it2), 0.0,                 // CVT_PT2( __LL3_it2)
	CVT_PT2( __LL3_it2), 0.0,                 // CVT_PT2(-__LL3_it2)
	CVT_PT2( __LL3_it1), 0.0,                 // CVT_PT2(-__LL3_it1)
};

static double QE1_it [] = {
	-__LL1_it, -1,         __LL1_it, -1,
	 1,        -__LL1_it,  1,         __LL1_it,
	-__LL1_it,  1,         __LL1_it,  1,
	-1,        -__LL1_it, -1,         __LL1_it,
};
static double QE2_it [] = {
	-__LL2_it, -1,         0, -1,  __LL2_it, -1,
	 1,         __LL2_it,  1,  0,  1,         __LL2_it,
	-__LL2_it,  1,         0,  1,  __LL2_it,  1,
	-1,         __LL2_it, -1,  0, -1,         __LL2_it,
};
static double QE3_it [] = {
	-__LL3_it1, -1,         -__LL3_it2, -1,          __LL3_it2, -1,          __LL3_it1, -1,
	-1,         -__LL3_it1, -1,         -__LL3_it1, -1,          __LL3_it1, -1,          __LL3_it1,
	-__LL3_it1,  1,         -__LL3_it2,  1,          __LL3_it2,  1,          __LL3_it1,  1,
	 1,         -__LL3_it1,  1,         -__LL3_it1,  1,          __LL3_it1,  1,          __LL3_it1,
};

//MODULE SHAPE_FUNC
/*******************************************************************************
 *
 * �e��`��֐��̒�`
 *
 ******************************************************************************/
/* ���v�f�̌`��֐��B�͈͂�[-1,1]
 * ��ϐ��p�̊֐��̒�`���s���A���̏�ő��ϐ��p�̌`��֐��̃��b�p�[���`����B
 */
static void __LL1_shape_func (double *N, const double t)
{
	N [0] = (1-t)/2;
	N [1] = (1+t)/2;
}

static void __LL2_shape_func (double *N, const double t)
{
	N [0] = t*(t-1)/2;
	N [1] = (1-t*t);
	N [2] = t*(t+1)/2;
}

static void __LL3_shape_func (double *N, const double t)
{
	N [0] = (1-t)*(9*t*t-1)/16;
	N [1] = 9*(1-t*t)*(1-3*t)/16;
	N [2] = 9*(1-t*t)*(1+3*t)/16;
	N [3] = (1+t)*(9*t*t-1)/16;
}

static void __dLL1_shape_func (double *N, const double t)
{
	N [0] = -0.5;
	N [1] =  0.5;
}

static void __dLL2_shape_func (double *N, const double t)
{
	N [0] = (2*t-1)/2;
	N [1] = -2*t;
	N [2] = (2*t+1)/2;
}

static void __dLL3_shape_func (double *N, const double t)
{
	N [0] = -(27*t*t-18*t-1)/16;
	N [1] =  9*(9*t*t-2*t-3)/16;
	N [2] = 9*(-9*t*t-2*t+3)/16;
	N [3] =  (27*t*t+18*t-1)/16;
}

/* ���ϐ��p�̃��b�p�[�B t�Ŏg����͍̂ŏ��̍��ڂł���B*/
void LL1 (double *N, const double *t) { __LL1_shape_func(N,*t); }
void LL2 (double *N, const double *t) { __LL2_shape_func(N,*t); }
void LL3 (double *N, const double *t) { __LL3_shape_func(N,*t); }
void dLL1 (double *N, const double *t) { __dLL1_shape_func(N,*t); }
void dLL2 (double *N, const double *t) { __dLL2_shape_func(N,*t); }
void dLL3 (double *N, const double *t) { __dLL3_shape_func(N,*t); }

/* �O�p�`�v�f�̌`��֐� */
void Tr1 (double *N, const double *t) {
	N [0] = t [0];
	N [1] = t [1];
	N [2] = 1 - N [0] - N [1];
}

void dTr1 (double *N, const double *t) {
	N [0] =  1;
	N [1] =  0;
	N [2] = -1;
	N [3] =  0;
	N [4] =  1;
	N [5] = -1;
}

void Tr2 (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [0] = t1*(2*t1-1);
	N [1] = t2*(2*t2-1);
	N [2] = t3*(2*t3-1);
	N [3] = 4*t1*t2;
	N [4] = 4*t2*t3;
	N [5] = 4*t3*t1;
}

void dTr2 (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [ 0] = 4*t1-1;
	N [ 1] = 0;
	N [ 2] = -4*t3+1;
	N [ 3] = 4*t2;
	N [ 4] = -4*t2;
	N [ 5] = 4*(t3-t1);

	N [ 6] = 0;
	N [ 7] = 4*t2-1;
	N [ 8] = -4*t3+1;
	N [ 9] = 4*t1;
	N [10] = 4*(t3-t2);
	N [11] = -4*t1;
}

void Tr3L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [0] = t1 * (3 * t1 - 2) * (3 * t1 - 1) / 2; // 1,0,0
	N [1] = t2 * (3 * t2 - 2) * (3 * t2 - 1) / 2; // 0,1,0
	N [2] = t3 * (3 * t3 - 2) * (3 * t3 - 1) / 2; // 0,0,1
	N [3] = 9 * t1 * t2 * (3 * t1 - 1) / 2;  // 2/3 1/3 0
	N [4] = 9 * t1 * t2 * (3 * t2 - 1) / 2;  // 1/3 2/3 0
	N [5] = 9 * t2 * t3 * (3 * t2 - 1) / 2;  // 0   2/3 1/3
	N [6] = 9 * t2 * t3 * (3 * t3 - 1) / 2;  // 0   1/3 2/3
	N [7] = 9 * t3 * t1 * (3 * t3 - 1) / 2;  // 1/3 0   2/3
	N [8] = 9 * t3 * t1 * (3 * t1 - 1) / 2;  // 2/3 0   1/3
	N [9] = 27 * t1 * t2 * t3;   // 1/3 1/3 1/3
}

void dTr3L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [ 0] =  (27 * t1 *t1 -18 * t1 + 2) / 2;
	N [ 1] =  0;
	N [ 2] = -(27 * t3 *t3 -18 * t3 + 2) / 2;
	N [ 3] =  9 * t2 * (6 * t1 - 1) / 2;
	N [ 4] =  9 * t2 * (3 * t2 - 1) / 2;
	N [ 5] = -9 * t2 * (3 * t2 - 1) / 2;
	N [ 6] = -9 * t2 * (6 * t3 - 1) / 2;
	N [ 7] = 9 * (3 * t3 * t3 - 6 * t1 * t3 - t3 + t1) / 2;
	N [ 8] = 9 * (6 * t1 * t3 - t3 - 3 * t1 * t1 + t1)/ 2;
	N [ 9] = 27 * t2 * (t3 - t1);

	N [10] = 0;
	N [11] = (27 * t2 *t2 -18 * t2 + 2) / 2;
	N [12] = -(27 * t3 *t3 -18 * t3 + 2) / 2;
	N [13] = 9 * t1 * (3 * t1 - 1) / 2;
	N [14] = 9 * t1 * (6 * t2 - 1) / 2;
	N [15] = 9 * (-3 * t2 * t2 + t2 + 6 * t2 * t3 - t3) / 2;
	N [16] = 9 * (3 * t3 * t3 - t3 -6 * t2 * t3 + t2) / 2;
	N [17] = -9 * t1 * (6 * t3 - 1) / 2;
	N [18] = -9 * t1 * (3 * t1 - 1) / 2;
	N [19] = 27 * t1 * (t3 - t2);
}

/* �����ɐߓ_���Ȃ�3���O�p�`�v�f�̌`��֐��B */
void Tr3S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [0] = t1 * (3 * t1 - 2) * (3 * t1 - 1) / 2;
	N [1] = t2 * (3 * t2 - 2) * (3 * t2 - 1) / 2;
	N [2] = t3 * (3 * t3 - 2) * (3 * t3 - 1) / 2;
	N [3] = 9 * t1 * t2 * (2 * t1 - t2) / 2;
	N [4] = 9 * t1 * t2 * (2 * t2 - t1) / 2;
	N [5] = 9 * t2 * t3 * (2 * t2 - t3) / 2;
	N [6] = 9 * t2 * t3 * (2 * t3 - t2) / 2;
	N [7] = 9 * t3 * t1 * (2 * t3 - t1) / 2;
	N [8] = 9 * t3 * t1 * (2 * t1 - t3) / 2;
}

void dTr3S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1], t3 = 1-t1-t2;
	N [ 0] =  (27 * t1 *t1 -18 * t1 + 2) / 2;
	N [ 1] =  0;
	N [ 2] = -(27 * t3 *t3 -18 * t3 + 2) / 2;
	N [ 3] = 9 * t2 * (4 * t1 - t2) / 2;
	N [ 4] = 9 * t2 * (2 * t2 - 2 * t1) / 2;
	N [ 5] = 9 * t2 * (-2 * t2 + 2 * t3) / 2;
	N [ 6] = 9 * t2 * (-4 * t3 + t2) / 2;
	N [ 7] = 9 * (2 * t3 * t3 - 6 * t1 * t3 + t1 * t1) / 2;
	N [ 8] = 9 * (6 * t3 * t1 - t3 * t3 - 2 * t1 * t1) / 2;

	N [ 9] = 0;
	N [10] = (27 * t2 *t2 -18 * t2 + 2) / 2;
	N [11] = -(27 * t3 *t3 -18 * t3 + 2) / 2;
	N [12] = 9 * t1 * (2 * t1 - 2 * t2) / 2;
	N [13] = 9 * t1 * (4 * t2 - t1) / 2;
	N [14] = 9 * (6 * t3 * t2 - t3 * t3 - 2 * t2 * t2) / 2;
	N [15] = 9 * (2 * t3 * t3 - 6 * t2 * t3 + t2 * t2) / 2;
	N [16] = 9 * t1 * (-4 * t3 + t1) / 2;
	N [17] = 9 * t1 * (-2 * t1 + 2 * t3) / 2;
}

/* �l�p�`�v�f�̌`��֐��B
 * 2���ȏ�̃Z�����f�B�s�e�B�v�f�̒�`�̓c�B�G���L�[���B�b�c�̒�{�ɋ������B
 * ���O�����W�F�v�f�͐��v�f�̌`��֐����|�����킹��Ƃł���B
 */
static void __Q1 (double *N, const double *L1, const double *L2) {
	N [0] = L1 [0] * L2 [0];
	N [1] = L1 [1] * L2 [0];
	N [2] = L1 [1] * L2 [1];
	N [3] = L1 [0] * L2 [1];
}

static void __Q2L (double *N, const double *L1, const double *L2) {
	N [0] = L1 [0] * L2 [0];
	N [1] = L1 [2] * L2 [0];
	N [2] = L1 [2] * L2 [2];
	N [3] = L1 [0] * L2 [2];

	N [4] = L1 [1] * L2 [0];
	N [5] = L1 [2] * L2 [1];
	N [6] = L1 [1] * L2 [2];
	N [7] = L1 [0] * L2 [1];

	N [8] = L1 [1] * L2 [1];
}

static void __Q2S (double *N, const double *L1, const double *L2, const double *L21, const double *L22) {
	N [4] = L21 [1] * L2 [0];
	N [5] = L1 [1] * L22 [1];
	N [6] = L21 [1] * L2 [1];
	N [7] = L1 [0] * L22 [1];
	N [0] = L1 [0] * L2 [0] - (N [4] + N [7]) / 2;
	N [1] = L1 [1] * L2 [0] - (N [5] + N [4]) / 2;
	N [2] = L1 [1] * L2 [1] - (N [6] + N [5]) / 2;
	N [3] = L1 [0] * L2 [1] - (N [7] + N [6]) / 2;
}

static void __Q3L (double *N, const double *L1, const double *L2) {
	N [ 0] = L1 [0] * L2 [0];
	N [ 1] = L1 [3] * L2 [0];
	N [ 2] = L1 [3] * L2 [3];
	N [ 3] = L1 [0] * L2 [3];

	N [ 4] = L1 [1] * L2 [0];
	N [ 5] = L1 [2] * L2 [0];
	N [ 6] = L1 [3] * L2 [1];
	N [ 7] = L1 [3] * L2 [2];
	N [ 8] = L1 [2] * L2 [3];
	N [ 9] = L1 [1] * L2 [3];
	N [10] = L1 [0] * L2 [2];
	N [11] = L1 [0] * L2 [1];

	N [12] = L1 [1] * L2 [1];
	N [13] = L1 [2] * L2 [1];
	N [14] = L1 [2] * L2 [2];
	N [15] = L1 [1] * L2 [2];
}

static void __Q3S (double *N, const double *L1, const double *L2, const double *L31, const double *L32) {
	N [ 4] = L31 [1] * L2  [0];
	N [ 5] = L31 [2] * L2  [0];
	N [ 6] = L1  [1] * L32 [1];
	N [ 7] = L1  [1] * L32 [2];
	N [ 8] = L31 [2] * L2  [1];
	N [ 9] = L31 [1] * L2  [1];
	N [10] = L1  [0] * L32 [2];
	N [11] = L1  [0] * L32 [1];
	N [0] = L1 [0] * L2 [0] - 2 * (N [11] + N [ 4]) / 3 - (N [10] + N [ 5]) / 3;
	N [1] = L1 [1] * L2 [0] - 2 * (N [ 5] + N [ 6]) / 3 - (N [ 4] + N [ 7]) / 3;
	N [2] = L1 [1] * L2 [1] - 2 * (N [ 7] + N [ 8]) / 3 - (N [ 6] + N [ 9]) / 3;
	N [3] = L1 [0] * L2 [1] - 2 * (N [ 9] + N [10]) / 3 - (N [ 8] + N [11]) / 3;
}

void Q1 (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__Q1 (N, L1, L2);
}

void dQ1 (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	double dL1 [2], dL2 [2];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__dLL1_shape_func (dL1, t1);
	__dLL1_shape_func (dL2, t2);
	__Q1 (N + 0, dL1, L2);
	__Q1 (N + 4, L1, dL2);
}

void Q2L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [3], L2 [3];
	__LL2_shape_func (L1, t1);
	__LL2_shape_func (L2, t2);
	__Q2L (N, L1, L2);
}

void dQ2L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [3], L2 [3];
	double dL1 [3], dL2 [3];
	__LL2_shape_func (L1, t1);
	__LL2_shape_func (L2, t2);
	__dLL2_shape_func (dL1, t1);
	__dLL2_shape_func (dL2, t2);
	__Q2L (N + 0, dL1, L2);
	__Q2L (N + 9, L1, dL2);
}

void Q2S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	double L21 [3], L22 [3];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__LL2_shape_func (L21, t1);
	__LL2_shape_func (L22, t2);
	__Q2S (N, L1, L2, L21, L22);
}

void dQ2S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	double L21 [3], L22 [3];
	double dL1 [2], dL2 [2];
	double dL21 [3], dL22 [3];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__LL2_shape_func (L21, t1);
	__LL2_shape_func (L22, t2);
	__dLL1_shape_func (dL1, t1);
	__dLL1_shape_func (dL2, t2);
	__dLL2_shape_func (dL21, t1);
	__dLL2_shape_func (dL22, t2);
	__Q2S (N + 0, dL1, L2, dL21, L22);
	__Q2S (N + 8, L1, dL2, L21, dL22);
}

void Q3L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [4], L2 [4];
	__LL3_shape_func (L1, t1);
	__LL3_shape_func (L2, t2);
	__Q3L (N, L1, L2);
}

void dQ3L (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [4], L2 [4];
	double dL1 [4], dL2 [4];
	__LL3_shape_func (L1, t1);
	__LL3_shape_func (L2, t2);
	__dLL3_shape_func (dL1, t1);
	__dLL3_shape_func (dL2, t2);
	__Q3L (N +  0, dL1, L2);
	__Q3L (N + 16, L1, dL2);
}

void Q3S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	double L31 [4], L32 [4];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__LL3_shape_func (L31, t1);
	__LL3_shape_func (L32, t2);
	__Q3S (N, L1, L2, L31, L32);
}

void dQ3S (double *N, const double *t) {
	double t1 = t[0], t2 = t[1];
	double L1 [2], L2 [2];
	double L31 [4], L32 [4];
	double dL1 [2], dL2 [2];
	double dL31 [4], dL32 [4];
	__LL1_shape_func (L1, t1);
	__LL1_shape_func (L2, t2);
	__LL3_shape_func (L31, t1);
	__LL3_shape_func (L32, t2);
	__dLL1_shape_func (dL1, t1);
	__dLL1_shape_func (dL2, t2);
	__dLL3_shape_func (dL31, t1);
	__dLL3_shape_func (dL32, t2);
	__Q3S (N +  0, dL1, L2, dL31, L32);
	__Q3S (N + 12, L1, dL2, L31, dL32);
}

/*
 * �ʏ��f�[�^�[
 */
// static int _FINFO_L1 [] = { 0, 1, };
// static int _FINFO_L2 [] = { 0, 1, 2, };
// static int _FINFO_L3 [] = { 0, 1, 2, 3, };
static int _FINFO_Tr1 [] = {
	0, 1,  1, 2,  2, 0,
};
static int _FINFO_Tr2 [] = {
	0, 3, 1,  1, 4, 2,  2, 5, 0,
};
static int _FINFO_Tr3 [] = {
	0, 3, 4, 1,  1, 5, 6, 2,  2, 0, 7, 8, 0,
};
static int _FINFO_Q1 [] = {
	0, 1,  1, 2,  2, 3,  3, 0,
};
static int _FINFO_Q2 [] = {
	0, 4, 1,  1, 5, 2,  2, 6, 3,  3, 7, 0,
};
static int _FINFO_Q3 [] = {
	0, 4, 5, 1,  1, 6, 7, 2,  2, 8, 9, 3,  3, 10,11, 0,
};

//

//MODULE DATA_DEF
/*******************************************************************************
 *
 * �L���v�f��͂ɕK�v�ȃf�[�^�[�^�̒�`
 *
 ******************************************************************************/


/* �����l�f�[�^�[�B��{�^�̂ق��A�����̍\����͂̂��߂̃f�[�^�[���`����B */
typedef struct MatData *pMatData;
typedef struct {
	int number;
	char *name;
	int type;
	pMatData ref_mat;	// �Q�ƕ����e�[�u��
	double init_temp;	// �������x
	double density;		// ���ʖ��x
} MatData;

typedef struct {
	int number;
	char *name;
	int type;
	pMatData ref_mat;
	double init_temp;
	double density;
	double young, poisson, shar_mod;	// �����O���A�|�A�\����A���c�e���W��
} IsoMatData;

/* �􉽓����l */
typedef struct {
	int number;
	char *name;
	int type;
} GeomData;

typedef GeomData PSN_GeomData;
typedef GeomData AXSol_GeomData;

typedef struct {
	int number;
	char *name;
	int type;
	double thick;
} PSS_GeomData;

/* �p�[�c�f�[�^�[�B�����l�Ɗ􉽓����l�����B�v�f����Q�Ƃ����B */
typedef struct {
	int number;
	char *name;
	MatData *mat;
	GeomData *geom;
	double *init_mat; // ���������}�g���b�N�X
} PartData;

/* �ߓ_�f�[�^�[ */
typedef struct {
	int number;
	double *coord;
} Node;

/* �v�f�^�C�v�\���� */
typedef struct {
	int number;
	char *name;
	int prob, ndim, ndof, nvstat, ijh_scale, ijh_add;
} ElementTypeInfo;

/* �`��֐��̃|�C���^ */
typedef void (*SF)(double *, const double *);

/* �v�f�ϕ��v�Z���\���̂̒�` */
struct ElementInfo;
typedef struct {
	int number;
	int nnode, norder, ntint, ndim; // �v�f������̐ߓ_���A�����A�S�ϕ��_���A������
	SF sf, dsf;   // �`��֐��ւ̃|�C���^�[
	double *N, *dN, *ipt, *wt; // �`��֐��ƌ`��֐��̔����`�A�ϕ��_���Əd�݂̔z��
} ElemIntegralInfo;

/* �v�f���\���̂̒�`�B�e�v�f����Q�Ƃ����B
 * �v�Z��񂪓���邪�A�v�f�S�̂Ɗe�\�ʂ̌`��֐����܂ނ��߂ł���B
 * �v�f�\�ʂ̌v�Z�͕��z�׏d�̌v�Z�ɗp����B
 */
typedef struct {
	int number;
	char *name;
	int nvert, nface;  // �v�f������̐ߓ_���A���_(���ߓ_)�̐��A�\�ʂ̐�
	ElemIntegralInfo *info1, *info2, *info3, *info4; // �v�f�S�̂ƕ\�ʂ̌v�Z���
	int *edge_list;
} ElementInfo;

/* �v�Z���f�[�^�[�B�\����͂ɓ��������f�[�^�[�̑g�ݕ������Ă��邪�A����ȊO�ɂ��g�p�\�ł���B */
typedef struct {
	double detJ, sigy, *iJH;  // ���R�r�A���̍s��l��J-1*H�̌��ʁB
	double *stress, *strain, *FT; // ���͂ƂЂ��݁A�ό`���z
	double *stat;   // ��ԕϐ�(���x�Ȃ�)
} ElementCalcInfo;

/* �v�f�f�[�^�[ */
typedef struct {
	int number;
	PartData *part;
	ElementInfo *info;
	int *conn;
	ElementCalcInfo *cinfo;
} Element;

/* �Z�b�g�f�[�^�[ */
/* �ŏ��ɃZ�b�g�Ɋ܂܂��f�[�^�[�^���`����B
 * ���̂悤�Ȃ��̂͐��l���̂��̂��g�p�����A�V���{���萔���g�p����B
 */
enum { SET_NONE=0, SET_NODE, SET_SEG, SET_ELEM, SET_FACE, SET_PART, };
typedef struct {
	int number;
	char *name;
	int type;
	int ndata, *data;
} SetData;

/*
 * ���E�����̒�`
 */

/* �S�������̒�` */
typedef struct {
	int number;
	char *name;
	int set_type;		// �Z�b�g�^�C�v�BSET_NODE�ŌŒ�
	int ndata, *data;
	int flags [6];		// �S��������������ꍇ�A1���Z�b�g����B
	double value [6];	// �S�������̒l�B
} FixDisp;

/* �׏d�����̒�` */
// �׏d��������ΏہB�ߓ_�E�\�ʁE�̐�(�v�f�S��)
enum { BC_POINT=1, BC_BODY = SET_ELEM, BC_FACE = SET_FACE, };

// ���͕��@�B
// ����(0)�A�e�����̒l�A�����Ƒ傫��(�����͐��K�������)�A�����Ƒ傫���̐ρA������������(�ӁE�ʉ׏d)
enum { BC_NONE, BC_VAL, BC_VEC, BC_VAL_VEC, BC_NORMAL, };

// �׏d�����̒ǉ��f�[�^�[�B���͎���10�{����
enum { BC_LENGTH, BC_AREA = 1, };		// �P�ʒ���(�����̂�)�ƒP�ʖʐς��Ƃ̓��́B
						// �����ɉ׏d��������v�f�ŗL���B�\���b�h�V�F���E�ő̗v�f�ł͖��������B
enum { BC_VOLVAL, BC_DENSITY = 1, };		// �̐ϗ͂̃t���O�B�f�t�H���g�͔�d���l��

typedef struct {
	int number;
	char *name;
	int set_type;		// �Z�b�g�^�C�v�B
	int ndata, *data;
	int load_type;
	int val_type;	// ���3�̗񋓑̗v�f�̊i�[��
	double value [4], load [3];	// ���͒l�Ǝ��ۂ̉׏d
} LoadData;

//MODULE GRLBAL_DATA
/*****************************************************************************
 * ���ϐ��̐錾
 *****************************************************************************/

/* ����Ԃł̍��W�n�̐��Ǝ��R�x�̐� */
int mdof = 2;
int mdim = 2;

/* ��̓^�C�g�� */
char *title;

/* �ėp���e�l */
double tolerance = 1e-10;

/* ���ԍ��B�񋓑̂Ƃ��邱�Ƃŕ�����₷�������B */
enum { PT_None=0, PT_PSS, PT_PSN, PT_AXSOL, };
int prob;

/* �p�[�c�╨���l�A�􉽓����l�̔z�� */
int npart, nmat, ngeom;
PartData *part;
MatData **mat;
GeomData **geom;

/* �ߓ_�Ɨv�f�̔z��̐錾 */
int ntnode, ntelem;
Node *node;
Element *elem;
double *coord;

/* �Z�b�g�f�[�^�[�̐錾 */
int nset;
SetData *set;

/* ���E�����f�[�^�[�̐錾 */
int nfdisp, npload, nfload, nbload;
FixDisp *fdisp;
LoadData *pload, *fload, *bload;

/* �S�̍����}�g���b�N�X�̐��l�i�[���@�F�����A�t���}�g���b�N�X�A�o���h�}�g���b�N�X�A�Ώ̍s��(�㑤�̂�) */
enum { RANGE_NONE, RANGE_ALL, RANGE_BEGIN_END, RANGE_SYMM_UPPER, };

/* �S�̍����}�g���b�N�X�ւ̐��l�̊i�[���@�A��L�̒l�̂����ꂩ������B */
int sys_range_mode;

/* �S���̏�ԁB���ɁA�S���Ȃ��A0�ȊO�ōS���A0�ōS���A�\�����Ɛڑ����Ă��Ȃ��ߓ_ */
enum { KLDOF_NONE, KLDOF_NONZERO, KLDOF_ZERO, KLDOF_ORPHAN, };

/* ���C�����\���́B�S�̍����̏����i�[���� */
typedef struct {
	int pos, start, end, size;	/* �n�_�ƏI�_�̎��̈ʒu�ƃT�C�Y */
	double *p;
} K_line_info;
K_line_info *line_info;

/* �S�̍����̎��f�[�^�[�ƕψʂƉ׏d */
double *sysk, *lhs_value, *rhs_value;

/* �v�Z���ʂ��i�[ */
double *disp, *force;

/* �ߓ_�ԍ��Ǝ��R�x����S�̍����̃C���f�b�N�X���t�Ɉ������X�g */
int *rev_line_info, *decode_rev_info;
/* �e���E�����̊J�n�ʒu�ƏI���ʒu */
int rank_line_info [4][2];

/* �ߓ_���Q�Ƃ���v�f���Ɗe���R�x�̃^�C�v */
static int *nodal_dof_list;

//MODULE INIT_ELEMENT_INFO
/*
 * �v�f���e�[�u���̒�`
*/
ElementTypeInfo etinfo [] = {
	/* number   name                prob    dim dof stat scl add */
	{ 0,        "None",             PT_None,  0,  0,   0,  0, 0, },
	{ PT_PSS,   "Plane stress",     PT_PSS,   2,  2,   3,  2, 1, },
	{ PT_PSN,   "Plane strain",     PT_PSN,   2,  2,   3,  2, 0, },
	{ PT_AXSOL, "Aximentric solid", PT_AXSOL, 2,  2,   4,  3, 1, },
};

ElemIntegralInfo eiinfo [] = {
/*        num nod ord int dim sf1   sf2    N1    dN1   ipt1    wt1     */
	{  0,  0, 0,  0,  0,  NULL, NULL,  NULL, NULL, NULL,   NULL,   },
	{  1,  2, 1,  2,  1,  LL1,  dLL1,  NULL, NULL, LL1_it, LL1_wt, },
	{  2,  3, 2,  3,  1,  LL2,  dLL2,  NULL, NULL, LL2_it, LL2_wt, },
	{  3,  4, 3,  4,  1,  LL3,  dLL3,  NULL, NULL, LL3_it, LL3_wt, },
	{  4,  3, 1,  1,  2,  Tr1,  dTr1,  NULL, NULL, Tr1_it, Tr1_wt, },
	{  5,  6, 2,  3,  2,  Tr2,  dTr2,  NULL, NULL, Tr2_it, Tr2_wt, },
	{  6, 10, 3,  7,  2,  Tr3L, dTr3L, NULL, NULL, Tr3_it, Tr3_wt, },
	{  7,  9, 3,  7,  2,  Tr3S, dTr3S, NULL, NULL, Tr3_it, Tr3_wt, },
	{  8,  4, 1,  4,  2,  Q1,   dQ1,   NULL, NULL, Q1_it,  Q1_wt,  },
	{  9,  9, 2,  9,  2,  Q2L,  dQ2L,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 10,  8, 2,  9,  2,  Q2S,  dQ2S,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 11, 16, 3, 16,  2,  Q3L,  dQ3L,  NULL, NULL, Q3_it,  Q3_wt,  },
	{ 12, 12, 3, 16,  2,  Q3S,  dQ3S,  NULL, NULL, Q3_it,  Q3_wt,  },
	{ 13,  6, 2,  1,  2,  Tr2,  dTr2,  NULL, NULL, Tr1_it, Tr1_wt, },
	{ 14, 10, 3,  3,  2,  Tr3L, dTr3L, NULL, NULL, Tr2_it, Tr2_wt, },
	{ 15,  9, 3,  3,  2,  Tr3S, dTr3S, NULL, NULL, Tr2_it, Tr2_wt, },
	{ 16,  4, 1,  1,  2,  Q1,   dQ1,   NULL, NULL, Q0_it,  Q0_wt,  },
	{ 17,  9, 2,  4,  2,  Q2L,  dQ2L,  NULL, NULL, Q1_it,  Q1_wt,  },
	{ 18,  8, 2,  4,  2,  Q2S,  dQ2S,  NULL, NULL, Q1_it,  Q1_wt,  },
	{ 19, 16, 3,  9,  2,  Q3L,  dQ3L,  NULL, NULL, Q2_it,  Q2_wt,  },
	{ 20, 12, 3,  9,  2,  Q3S,  dQ3S,  NULL, NULL, Q2_it,  Q2_wt,  },
};

ElementInfo einfo [] = {
/*        num name    vert face info1         info2  info3        info4         edge       */
	{  0, NULL,      0,  0,  NULL,         NULL, NULL,       NULL,          NULL,       },
	{  1, "Tria1",   3,  3,  &eiinfo [ 4], NULL, &eiinfo [1], &eiinfo [ 4], _FINFO_Tr1, },
	{  2, "Tria2",   3,  3,  &eiinfo [ 5], NULL, &eiinfo [2], &eiinfo [ 5], _FINFO_Tr2, },
	{  3, "Tria3L",  3,  3,  &eiinfo [ 6], NULL, &eiinfo [3], &eiinfo [ 6], _FINFO_Tr3, },
	{  4, "Tria3S",  3,  3,  &eiinfo [ 7], NULL, &eiinfo [3], &eiinfo [ 7], _FINFO_Tr3, },
	{  5, "Quad1",   4,  4,  &eiinfo [ 8], NULL, &eiinfo [1], &eiinfo [ 8], _FINFO_Q1,  },
	{  6, "Quad2L",  4,  4,  &eiinfo [ 9], NULL, &eiinfo [2], &eiinfo [ 9], _FINFO_Q2,  },
	{  7, "Quad2S",  4,  4,  &eiinfo [10], NULL, &eiinfo [2], &eiinfo [10], _FINFO_Q2,  },
	{  8, "Quad3L",  4,  4,  &eiinfo [11], NULL, &eiinfo [3], &eiinfo [11], _FINFO_Q3,  },
	{  9, "Quad3S",  4,  4,  &eiinfo [12], NULL, &eiinfo [3], &eiinfo [12], _FINFO_Q3,  },
	{ 10, "Tria2R",  3,  3,  &eiinfo [13], NULL, &eiinfo [2], &eiinfo [ 5], _FINFO_Tr2, },
	{ 11, "Tria3LR", 3,  3,  &eiinfo [14], NULL, &eiinfo [3], &eiinfo [ 6], _FINFO_Tr3, },
	{ 12, "Tria3SR", 3,  3,  &eiinfo [15], NULL, &eiinfo [3], &eiinfo [ 7], _FINFO_Tr3, },
	{ 13, "Quad1R",  4,  4,  &eiinfo [16], NULL, &eiinfo [1], &eiinfo [ 8], _FINFO_Q1,  },
	{ 14, "Quad2LR", 4,  4,  &eiinfo [17], NULL, &eiinfo [2], &eiinfo [ 9], _FINFO_Q2,  },
	{ 15, "Quad2SR", 4,  4,  &eiinfo [18], NULL, &eiinfo [2], &eiinfo [10], _FINFO_Q2,  },
	{ 16, "Quad3LR", 4,  4,  &eiinfo [19], NULL, &eiinfo [3], &eiinfo [11], _FINFO_Q3,  },
	{ 17, "Quad3SR", 4,  4,  &eiinfo [20], NULL, &eiinfo [3], &eiinfo [12], _FINFO_Q3,  },
};

int netinfo = sizeof (etinfo) / sizeof (etinfo [0]);
int neiinfo = sizeof (eiinfo) / sizeof (eiinfo [0]);
int neinfo  = sizeof (einfo)  / sizeof (einfo [0]);

/*
*/
#include "debug.c"

/*
 * �֐��̒�`
 */

/* �v�f�ϕ��v�Z���\���̂̏����� */
void init_einfo ()
{
		// 0�Ԗڂ̃G���g���[�͎g�p���Ȃ��̂ŁA1�𑫂���1�Ԗڂ���A�N�Z�X����悤�ɂ��Ă���B
	ElemIntegralInfo *eip = eiinfo + 1;
	int i, j;
	for (i = 1; i < neiinfo; i++, eip++) {
		int nnode = eip->nnode;
		int ndim = eip->ndim;
		int ntint = eip->ntint;
		double *ipt = eip->ipt;
		// �}�g���b�N�X�̈�̊��蓖��
		eip->N  = malloc (sizeof (double) * ntint * nnode);
		eip->dN = malloc (sizeof (double) * ntint * nnode * ndim);
		for (j = 0; j < eip->ntint; j++) {
			(*eip->sf)  (eip->N  + j * nnode,        ipt + j * ndim);
			(*eip->dsf) (eip->dN + j * nnode * ndim, ipt + j * ndim);
		}
	}
}

//MODULE MISC_FUNCTIONS
/*
 * �G���֐��Q
 */
/* �ԍ�����z���̈ʒu��T���֐��Q */
int search_einfo (int n)
{
	int i;
	for (i = 0; i < neinfo; i++)
		if (einfo [i].number == n) break;
	return i;
}

int search_part (int n)
{
	int i;
	for (i = 0; i < npart; i++)
		if (part [i].number == n) break;
	return i;
}

int search_node (int n)
{
	int i;
	for (i = 0; i < ntnode; i++)
		if (node [i].number == n) break;
	return i;
}

int search_elem (int n)
{
	int i;
	for (i = 0; i < ntelem; i++)
		if (elem [i].number == n) break;
	return i;
}

int search_mat (int n)
{
	int i;
	for (i = 0; i < nmat; i++)
		if (mat [i]->number == i) break;
	return i;
}

// MODULE IO_DATA
/*
 * IO���[�`���Q
 */
static int nline = 0;
static char buffer [201];
void free_mesh_data ();
void free_bc_data ();
void free_data ();
void free_elem_data (int index);

static void error (int code, char *msg)
{
	fprintf (stderr, "%s at %d\n", msg, nline);
	free_data ();
	exit (code);
}

static void read_line (FILE *fin)
{
	char *p;
	fgets (buffer, sizeof (buffer), fin), nline++;
	while (*buffer == '#')
		fgets (buffer, sizeof (buffer), fin), nline++;
	for (p = buffer; !*p && *p != '#'; p++)
		;
	if (*p == '#')
		for (; !*p; p++) *p = 0;
	for (p = buffer; *p; p++)
		if (*p < 0x20) *p = 0;
}

static void read_str (FILE *fin, char **s)
{
	read_line (fin);
	*s = malloc (strlen(buffer)+1);
	strcpy (*s, buffer);
}

/* ���X�g�ɓ��Y�f�[�^�[��ǉ� */
void attach_node (int *ip, int n)
{
	int ptr = search_node (n);
	if (ptr == ntnode) error (-99, "Illegal Node number");
	*ip = ptr;
}

void attach_elem (int *ip, int n)
{
	int ptr = search_elem (n);
	if (ptr == ntelem) error (-99, "Illegal elem number");
	*ip = ptr;
}

void attach_part (int *ip, int n)
{
	int ptr = search_part (n);
	if (ptr == npart) error (-99, "Illegal part number");
	*ip = ptr;
}

/**/
static void read_header (FILE *fin)
{
	read_str (fin, &title);
	read_line (fin);
	sscanf (buffer, "%d%d%d%d%d%d%d", &prob, &npart, &nmat, &ngeom, &ntnode, &ntelem, &nset);
	if (prob <= 0 && prob > PT_AXSOL) error (-2, "Illegal progrem setting.");
	if (npart < 0) error (-2, "Illegal part count.");
	if (nmat  < 0) error (-2, "Illegal material table count.");
	if (ngeom < 0) error (-2, "Illegal geometric table count.");
	if (ntnode < 0) error (-2, "Illegal node count");
	if (ntelem < 0) error (-2, "Illegal element count.");
	if (nset  < 0) error (-2, "Illegal set list count.");
}

static void write_header (FILE *fout)
{
	fprintf (fout, "\n\nEFEM/C ver1.60\n");
	fprintf (fout, "\tby H.Hoshino\n\n");
	fprintf (fout, "Analysis title\n");
	fprintf (fout, "\t%s\n\n", title);
	fprintf (fout, "Analysis type\n");
	fprintf (fout, "\t%s\n\n", etinfo [prob].name);
}

static void alloc_fe_data ()
{
	part = calloc (sizeof (PartData), npart);
	if (part == NULL) goto Error;
	mat  = calloc (sizeof (MatData *), nmat);
	if (mat == NULL)  goto Error;
	geom = calloc (sizeof (GeomData *), ngeom);
	if (geom == NULL) goto Error;
	node = calloc (sizeof (Node), ntnode);
	coord = node [0].coord = calloc (sizeof (double), ntnode * mdof);
	if (node == NULL) goto Error;
	elem = calloc (sizeof (Element), ntelem);
	if (elem == NULL) goto Error;
	set  = calloc (sizeof (SetData), nset);
	if (set == NULL) goto Error;
	return;
Error:
	free_mesh_data ();
	error  (-3, "No enough memory.");
}

static void read_mat (FILE *fin)
{
	int i;
	for (i = 0; i < nmat; i++) {
		MatData m;
		IsoMatData *mp;
		read_line (fin);
		sscanf (buffer, "%d%d%lg", &m.number, &m.type, &m.density);
		read_str (fin, &m.name);
		mp = malloc (sizeof (IsoMatData));
		*(MatData*)mp = m;
		read_line (fin);
		sscanf (buffer, "%lg%lg", &mp->young, &mp->poisson);
		mp->shar_mod = mp->young / 2 / (1 + mp->poisson);
		mat [i] = (MatData*)mp;
	}
}

static void write_mat (FILE *fout)
{
	int i;
	fprintf (fout, "Material table\n");
	fprintf (fout, "\tnumber of table: %d\n", nmat);
	for (i = 0; i < nmat; i++) {
		IsoMatData *p = (IsoMatData*)mat [i];
		fprintf (fout, "\n\t%-20s%d\n", "Table number: ", p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",           p->name);
		fprintf (fout, "\t\t%-20s%g\n", "Mass density:",   p->density);
		fprintf (fout, "\t\t%-20s%g\n", "Young modulas:",  p->young);
		fprintf (fout, "\t\t%-20s%g\n", "Poisson ratio:",  p->poisson);
	}
	fputc ('\n', fout);
}

static void read_geom (FILE *fin)
{
	int i;
	GeomData g, *gp;
	for (i = 0; i < nmat; i++) {
		read_line (fin);
		sscanf (buffer, "%d", &g.number);
		g.type = prob;
		read_str (fin, &g.name);
		if (prob == PT_PSS)
			gp = malloc (sizeof (PSS_GeomData));
		else
			gp = malloc (sizeof (GeomData));
		*gp = g;
		geom [i] = gp;
		if (prob == PT_PSS) {
			PSS_GeomData *pg = (PSS_GeomData *) gp;
			read_line (fin);
			sscanf (buffer, "%lg", &pg->thick);
		}
	}
}

static void write_geom (FILE *fout)
{
	int i;
	fprintf (fout, "Geometrical table\n");
	fprintf (fout, "\tnumber of table: %d\n", ngeom);
	for (i = 0; i < ngeom; i++) {
		GeomData *p = geom [i];
		fprintf (fout, "\n\t%-20s%d\n", "Table number:", p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",         p->name);
		fprintf (fout, "\t\t%-20s%s\n", "Type:",         etinfo [p->type].name);
		if (prob == PT_PSS) {
			PSS_GeomData *pg = (PSS_GeomData *) geom [i];
			fprintf (fout, "\t\t%-20s%g\n", "Thick:", pg->thick);
		}
	}
	fputc ('\n', fout);
}

static void read_part (FILE *fin)
{
	int i, j, m, g;
	PartData *p;
	for (i = 0; i < npart; i++) {
		p = part + i;
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &m, &g);
		read_str (fin, &p->name);
		for (j = 0; j < nmat; j++) {
			if (mat [j]->number == m) {
				p->mat = mat [i];
				break;
			}
		}
		if (j == nmat) error (-4, "Illegal material table number");
		for (j = 0; j < ngeom; j++) {
			if (geom [j]->number == m) {
				p->geom = geom [i];
				break;
			}
		}
		if (j == ngeom) error (-4, "Illegal geometrical table number");
	}
}

static void write_part (FILE *fout)
{
	int i;
	fprintf (fout, "Part information\n");
	fprintf (fout, "\tnumber of part: %d\n", npart);
	for (i = 0; i < npart; i++) {
		PartData *p = part + i;
		fprintf (fout, "\n\t%-20s%d\n", "Part Number ",       p->number);
		fprintf (fout, "\t\t%-20s%s\n", "Name:",              p->name);
		fprintf (fout, "\t\t%-20s%s\n", "Material table:",    p->mat->name);
		fprintf (fout, "\t\t%-20s%s\n", "Geometrical table:", p->geom->name);
	}
	fputc ('\n', fout);
}

static void read_node (FILE *fin)
{
	int i;
	double *dp;
	Node *p;
	for (i = 0, p = node, dp = coord; i < ntnode; i++, p++, dp += mdof) {
		p->coord = dp;
		read_line (fin);
		sscanf (buffer, "%d%lg%lg", &p->number, &p->coord [0], &p->coord[1]);
	}
}

static void write_node (FILE *fout)
{
	int i, j;
	double *dp;
	Node *p;
	char *sp [] = { "X-coord", "Y-coord", "Z-coord", };
	fprintf (fout, "%s\n" "\t%s : %d\n", "Node list", "number of node", ntnode);
	fprintf (fout, "%10s", "number");
	for (i = 0; i < mdim; i++)
		fprintf (fout, "%15s", sp [i]);
	for (i = 0, p = node; i < ntnode; i++, p++) {
		fprintf (fout, "%10d", p->number);
		for (j = 0, dp = p->coord; j < mdim; j++, dp++)
			fprintf (fout, "%15.6g", *dp);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
}

static int iJH_size (const Element *p)
{
	return etinfo [prob].ijh_scale * p->info->info1->nnode + etinfo [prob].ijh_add;
}

static void alloc_cinfo (Element *p, int ijh_size)
{
	int i;
	ElementCalcInfo *cp;
	int ntint = p->info->info1->ntint;
	p->cinfo = calloc (sizeof (ElementCalcInfo), ntint);
	if (p->cinfo == NULL) free_elem_data (p - elem);
	for (i = 0, cp = p->cinfo; i < ntint; i++, cp++) {
		cp->iJH = calloc (sizeof (double), ijh_size + 22);
		if (cp->iJH == NULL) free_elem_data (i);
		cp->stress = cp->iJH    + ijh_size;
		cp->strain = cp->stress + 6;
		cp->FT     = cp->strain + 6;
		cp->stat   = cp->FT     + 9;
	}
}

static void read_elem (FILE *fin)
{
	int i, j, ip, ie, ptr, pos;
	Element *p;
	char *cp;

	for (i = 0, p = elem; i < ntelem; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%d%n", &p->number, &ip, &ie, &pos);
		ptr = search_part (ip);
		if (ptr == npart) error (-5, "Illegal part table number");
		p->part = part + ptr;
		ptr = search_einfo (ie);
		if (ptr == neinfo) error (-5, "Illegal element infomation table number");
		p->info = einfo + ptr;
		p->conn = calloc (sizeof (int), p->info->info1->nnode);
		if (p->conn == NULL) free_elem_data (i);
		alloc_cinfo (p, iJH_size (p));
		cp = buffer + pos;
		for (j = 0; j < p->info->info1->nnode; j++) {
			int pos1, n;
			sscanf (cp, "%d%n", &n, &pos1);
			if (pos == pos1) {
				read_line (fin);
				cp = buffer;
				continue;
			}
			pos += pos1;
			cp += pos1;
			attach_node (p->conn + j, n);
		}
	}
}

static void write_elem (FILE *fout)
{
	int i, j, k, nnode, index;
	Element *p;
	fprintf (fout, "%s\n" "\t%s : %d\n\n", "Element list", "number of element", ntelem);
	fprintf (fout, "%10s%10s%10s %s\n", "number", "part", "type", "connectivity");
	for (i = 0, p = elem; i < ntelem; i++, p++) {
		fprintf (fout, "%10d%10d%10s", p->number, p->part->number, p->info->name);
		nnode = p->info->info1->nnode;
		if (nnode <= 8) {
			for (j = 0; j < nnode; j++)
				fprintf (fout, "%10d", node [p->conn [j]].number);
			fputc ('\n', fout);
		} else {
			for (j = index = 0; j < 8; j++, index++)
				fprintf (fout, "%10d", node [p->conn [index]].number);
			fputc ('\n', fout);
			nnode -= 8;
			for (j = 0; j < nnode / 10; j++) {
				fprintf (fout, "%10c", " ");
				for (k = 0; k < 10; k++, index++)
					fprintf (fout, "%10d", node [p->conn [index]].number);
				fputc ('\n', fout);
			}
			nnode %= 10;
			if (nnode > 0) {
				for (j = 0; j < nnode; j++, index++)
					fprintf (fout, "%10d", node [p->conn [index]].number);
				putc ('\n', fout);
			}
		}
	}
	fputc ('\n', fout);
}

static void read_int_list (FILE *fin, int ndata, int *ip)
{
	int pos, idata;
	char *cp;
	idata = 0;
	while (idata < ndata) {
		read_line (fin);
		cp = buffer;
		while (cp < buffer + strlen (buffer) && idata < ndata) {
			sscanf (cp, "%d%n", ip, &pos), cp += pos, idata++;
			//fprintf (stdout, "%5d%5d\n", *ip, pos);
			ip++;
		}
	}
}

void attach_entities (int type, int ndata, int *data)
{
	int i, n, *ip;
	switch (type) {
	case SET_NODE:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_node (ip, *ip);
		break;
	case SET_SEG:
		break;
	case SET_ELEM:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_elem (ip, *ip);
		break;
	case SET_FACE:
		for (i = 0, ip = data; i < ndata; i++, ip++) {
			attach_elem (ip, *ip);
			n = *ip++;
			if ((*ip <= 0) || (*ip > elem [n].info->nface))
				error (-6, "Illegal edge/face number");
			*ip -= 1;
		}
		break;
	case SET_PART:
		for (i = 0, ip = data; i < ndata; i++, ip++)
			attach_part (ip, *ip);
		break;
	default:
		error (-6, "Illegal set type");
	}
}

static void read_set (FILE *fin, int nset, SetData *p)
{
	int i, scale, *ip;
	for (i = 0; i < nset; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &p->type, &p->ndata);
		if ((p->type <= SET_NONE) || (p->type > SET_PART))
			error (-6, "Illegal set type");
		scale = 1;
		if (p->type == SET_FACE) scale = 2;
		read_str (fin, &p->name);
		p->data = calloc (sizeof (int), p->ndata * scale);
		read_int_list (fin, p->ndata * scale, p->data);
		attach_entities (p->type, p->ndata, p->data);
	}
}

static void write_part_list (FILE *fout, int ndata, const int *ip)
{
	int i;
	for (i = 0; i < ndata; i++) {
		fprintf (fout, "%10d", part [*ip++].number);
		if (i % 8 == 0) fputc ('\n', fout);
	}
	if (i % 8) fputs ("\n", fout);
}

static void write_node_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++)
			fprintf (fout, "%10d", node [*ip++].number);
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++)
			fprintf (fout, "%10d", node [*ip++].number);
		fputs ("\n", fout);
	}
}

static void write_elem_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++)
			fprintf (fout, "%10d", elem [*ip++].number);
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++)
			fprintf (fout, "%10d", elem [*ip++].number);
		fputs ("\n", fout);
	}
}

static void write_face_list (FILE *fout, int ndata, const int *ip)
{
	int i, j;
	const int ct = 8;
	for (i = 0; i < ndata / ct; i++) {
		for (j = 0; j < ct; j++) {
			fprintf (fout, "%10d:%d", elem [ip [0]].number, ip [1] + 1);
			ip += 2;
		}
		fputc ('\n', fout);
	}
	if (ndata % ct) {
		for (i *= ct; i < ndata; i++) {
			fprintf (fout, "%10d:%d", elem [ip [0]].number, ip [1] + 1);
			ip += 2;
		}
		fputc ('\n', fout);
	}
}

static void write_set (FILE *fout, int nset, SetData *p)
{
	int i;
	fprintf (fout, "%s\n",          "Set data");
	fprintf (fout, "\t%s : %d\n\n", "number of set", nset);
	for (i = 0; i < nset; i++, p++) {
		fprintf (fout, "\t%20s%d\n", "Table number:", p->number);
		fprintf (fout, "\t%20s%s\n", "Name:",         p->name);
		fprintf (fout, "\t%20s\n",   "Data list");
		switch (p->type) {
		case SET_NODE:
			fputs ("Node list", fout);
			write_node_list (fout, p->ndata, p->data);
			break;
		case SET_SEG:
			break;
		case SET_ELEM:
			fputs ("Element list", fout);
			write_elem_list (fout, p->ndata, p->data);
			break;
		case SET_FACE:
			fputs ("Boundary edge/face list", fout);
			write_face_list (fout, p->ndata, p->data);
			break;
		case SET_PART:
			fputs ("Part list", fout);
			write_part_list (fout, p->ndata, p->data);
			break;
		default:
			error (-6, "Illegal set type");
		}
	}
	fputc ('\n', fout);
}

void io_data (FILE *fin, FILE *fout)
{
	read_header   (fin);
	write_header  (fout);
	alloc_fe_data ();

	read_mat   (fin);
	read_geom  (fin);
	read_part  (fin);
	write_part (fout);
	write_mat  (fout);
	write_geom (fout);

	read_node  (fin);
	write_node (fout);
	read_elem  (fin);
	write_elem (fout);

	read_set  (fin, nset, set);
	write_set (fout, nset, set);
	fflush (fout);
}

/**/
static void read_bc_header (FILE *fin)
{
	read_line (fin);
	sscanf (buffer, "%d%d%d%d", &nfdisp, &npload, &nfload, &nbload);
	if (nfdisp < 0) error  (-7, "Illegal fixed displacement data count.");
	if (npload < 0) error  (-7, "Illegal point load data count.");
	if (nfload < 0) error  (-7, "Illegal face load data count.");
	if (nbload < 0) error  (-7, "Illegal body load data count.");
}

static void write_bc_header (FILE *fout)
{
	fprintf (fout, "%s\n\t%s: %d\n\t%s: %d\n\t%s: %d\n\t%s: %d\n\n",
		"Boundary condition tables",
		"Count of fixed disp table", nfdisp,
		"Count of point load table", npload,
		"Count of distributed load table", nfload,
		"Count of volume load table", nbload);
}

static void alloc_bc_data ()
{
	fdisp = calloc (sizeof (FixDisp), nfdisp);
	if (fdisp == NULL) goto Error;
	pload = calloc (sizeof (LoadData), npload);
	if (pload == NULL) goto Error;
	fload = calloc (sizeof (LoadData), nfload);
	if (fload == NULL) goto Error;
	bload = calloc (sizeof (LoadData), nbload);
	if (bload == NULL) goto Error;
	return;
Error:
	free_bc_data ();
	error (-7, "no enough memory");
}

static void read_fdisp (FILE *fin)
{
	int i, j, pos, *ip;
	char *cp;
	double *dp;
	FixDisp *p;

	for (i = 0, p = fdisp; i < nfdisp; i++, p++) {
		read_line (fin);
		sscanf (buffer, "%d%d%n", &p->number, &p->ndata, &pos);
		for (j = 0, ip = p->flags, cp = buffer + pos; j < etinfo [prob].ndof; j++, ip++) {
			sscanf (cp, "%d%n", ip, &pos);
			cp += pos;
		}
		read_str (fin, &p->name);
		read_line (fin);
		for (j = 0, dp = p->value, cp = buffer; j < etinfo [prob].ndof; j++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
		sscanf (buffer, "%lg%lg", &p->value [0], &p->value [1]);
		p->data = calloc (sizeof (int), p->ndata);
		read_int_list (fin, p->ndata, p->data);
		p->set_type = SET_NODE;
		for (j = 0; j < p->ndata; j++)
			attach_node (p->data + j, p->data [j]);
	}
}

static void set_load_value (FILE *fin, int nload, LoadData *p)
{
	char *cp;
	double *dp;
	int i, pos;

	read_line (fin);
	if ((p->val_type % 10== BC_VAL) || (p->val_type % 10 == BC_NORMAL && p->load_type % 10 == BC_FACE)) {
		for (i = 0, dp = p->value, cp = buffer; i < etinfo [prob].ndim; i++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
	} else if (p->val_type % 10 <= BC_VAL_VEC) {
		for (i = 0, dp = p->value, cp = buffer; i < etinfo [prob].ndim + 1; i++, dp++) {
			sscanf (cp, "%lg%n", dp, &pos);
			cp += pos;
		}
	} else {
		free_bc_data ();
		error (-8, "illegal point load type");
	}
}

static void read_load (FILE *fin, int nload, LoadData *p, int ltype)
{
	int i, j, ndata, scale, *ip;

	for (i = 0; i < nload; i++, p++) {
		memset (p, 0, sizeof (LoadData));
		p->load_type = ltype;
		read_line (fin);
		sscanf (buffer, "%d%d%d", &p->number, &ndata, &p->val_type);
		read_str (fin, &p->name);
		set_load_value (fin, i, p);
		p->ndata = ndata;
		p->data = calloc (sizeof (int), ndata);
		scale = 1;
		if (ltype == BC_FACE) scale = 2;
		read_int_list (fin, ndata * scale, p->data);
		attach_entities (ltype, p->ndata, p->data);
	}
}

static void read_bc (FILE *fin)
{
	read_fdisp (fin);
	read_load (fin, npload, pload, BC_POINT);
	read_load (fin, nfload, fload, BC_FACE);
	read_load (fin, nbload, bload, BC_BODY);
}

static void write_fdisp (FILE *fout)
{
	int i, j;
	FixDisp *p;
	char *sp[] = { "X-dir", "Y-dir", "Z-dir", "RX-dir", "RY-dir", "RZ-dir", };
	for (i = 0, p = fdisp; i < nfdisp; i++, p++) {
		fprintf (fout, "\n%20s%d\n%20s%s\n", "Table number:", p->number, "Name:", p->name);
		fprintf (fout, "%20s\n", "Constraint dof");
		for (j = 0; j <  etinfo [prob].ndof; j++)
			if (p->flags [j]) fprintf (fout, "%28s:%g\n", sp [j], p->value[j]);
		fprintf (fout, "Apply nodes\n");
		write_node_list (fout, p->ndata, p->data);
		fputc ('\n', fout);
	}
}

typedef void (*fp_write_int_list) (FILE *fout, int, const int *);

static void write_load (FILE *fout, int ndata, const LoadData *p, fp_write_int_list fp)
{
	int i, j;
	char *type1_sp  [] = { "None", "Direction", "Director and value", "Director mult value", };
	char *type2a_sp [] = { "None", "par length", "par area", };
	char *type2b_sp [] = { "None", "gravity", "uniform load", };
	char **type2_sp    = (p->set_type == BC_BODY)? type2b_sp: type2a_sp;
	char *type3_sp  [] = { "None", "at node", "at segment", "at face", "at body", "at part", };

	for (i = 0; i < ndata; i++, p++) {
		fprintf (fout, "%20s%d\n" "%20s%s\n", "Table number:", p->number, "Name:", p->name);
		fprintf (fout, "%20s:%s %s %s\n",
				"Load types", type1_sp [p->val_type % 10],
				type2_sp [p->val_type / 10], type3_sp [p->set_type]);
		fprintf (fout, "%20s", "Load value and vector:");
		for (j = 0; j <  etinfo [prob].ndof; j++)
			fprintf (fout, "%15g", p->value[j]);
		fputc ('\n', fout);
		fprintf (fout, "Apply entitis\n");
		(*fp) (fout, p->ndata, p->data);
		fputc ('\n', fout);
	}
}

static void write_bc (FILE *fout)
{
	fprintf (fout, "%s\n\t%s: %d\n\n", "Fixed displacement data", "Number of tables", nfdisp);
	write_fdisp (fout);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Point load data", "Number of tables", npload);
	write_load (fout, npload, pload, write_node_list);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Distributed load data", "Number of tables", nfload);
	write_load (fout, nfload, fload, write_face_list);
	fprintf (fout, "%s\n\t%s: %d\n\n", "Volume load data", "Number of tables", nbload);
	write_load (fout, nbload, bload, write_elem_list);
}

void io_bc (FILE *fin, FILE *fout)
{
	read_bc_header (fin);
	write_bc_header (fout);
	alloc_bc_data ();
	fflush (fout);
	read_bc (fin);
	write_bc (fout);
	fflush (fout);
}

//MODULE MAT
/*
 * �����s��̌v�Z
 */
void calc_De_Iso_PSS_Mat(double *De, const MatData *tab) {
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - p->poisson * p->poisson);
	double *D[] = { De, De+3, De+6, };
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = f;
	D[0][1] = D[1][0] = f * p->poisson;
	D[0][2] = D[1][2] = D[2][0] = D[2][1]= 0;
	D[2][2] = p->shar_mod;
}

void calc_De_Iso_PSN_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+3, De+6, };
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = f * (1 - p->poisson);
	D[0][1] = D[1][0] = f * p->poisson;
	D[0][2] = D[1][2] = D[2][0] = D[2][1]= 0;
	D[2][2] = p->shar_mod;
}

void calc_De_Iso_AXSol_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+4, De+8, De+12, };
	int i;
	for (i = 0; i < 16; i++) De [i] = 0;
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = D[2][2] = f * (1 - p->poisson);
	D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = f * p->poisson;
	D[3][3] = p->shar_mod;
}
void calc_De_Iso_Solid_Mat(double *De, const MatData *tab)
{
	IsoMatData *p = (IsoMatData *)tab;
	double f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	double *D[] = { De, De+6, De+12, De+18, De+24, De+30, };
	int i;
	for (i = 0; i < 36; i++) De [i] = 0;
	p->shar_mod=p->young / (1 + p->poisson) / 2;
	D[0][0] = D[1][1] = D[2][2] = f * (1 - p->poisson);
	D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = f * p->poisson;
	D[3][3] = D[4][4] = D[5][5] = p->shar_mod;
}

/* �e�p�[�c�̏����������v�Z���� */
typedef void (*DmatFP) (double *, const MatData *);
void calc_D_mat()
{
	int i, size = etinfo [prob].nvstat;
	DmatFP fp [] = { NULL, calc_De_Iso_PSS_Mat, calc_De_Iso_PSN_Mat, calc_De_Iso_AXSol_Mat, };
	size *= size;
	for (i = 0; i < npart; i++) {
		PartData *p = part + i;
		p->init_mat = calloc (sizeof (double), size);
		(*fp [prob]) (p->init_mat, p->mat);
	}
}

/* �e�Y����͂̂��߂̃��[�`���Q */
// �Ð��������߂�
// �Ð����͐������͂̕��ϒl�ɓ�����
void calc_deviatoric_stres (double *estres, double *HSS, const double *stres)
{
	int i;
	*HSS = (stres [0] + stres [1] + stres [2])/3;
	for (i = 0; i < 3; i++) {
		estres [i]     = stres [i] - *HSS;
		estres [i + 3] = stres [i + 3];
	}

}

// �s�ϗʂ����߂�
void calc_values (double *res, const double *s)
{
	res [0] = s [0] + s [1] + s [2];
	res [1] = s[0] * s[1] + s[1] * s[2] + s[2] * s[0] - s[3] * s[3] - s[4] * s[4] - s[5] * s[5];
	res [2] = s[0] * s[1] * s[2] + 2*s[3] * s[4] * s[5] - s[0] * s[3] * s[3] - s[1] * s[4] * s[4] - s[2] * s[5] * s[5];
}

double sign (double x)
{
	return x>=0? 1: -1;
}

// �剞�͂����߂�
void calc_principle_values (double *res, const double *s)
{
}

// �~�[�[�X�̑������͂����߂�
double mieses_stres(const double *s)
{
	return sqrt(s[0]*s[0]+s[1]*s[1]+s[2]*s[2]-s[0]*s[1]-s[1]*s[2]-s[2]*s[0]+3*(s[3]*s[3]+s[4]*s[4]+s[5]*s[5]));
}

// �Y�������}�g���b�N�X���v�Z����B
// �o�T�͎O�D[1976]�ɂ��
// ���ʉ��͂̂ݑ��ƈقȂ邪�A���͏�Ԃ̉ߒ����炱�̂悤�Ȍv�Z�ɂȂ��Ă���
void calc_PSS_IsoDp_Mat (double *_D, const MatData *tab, const double *stres, double hard)
{
	double *D[] = { _D, _D+3, _D+6, };
	double S[7], G, f, sb, sigbar, estres [6];
	IsoMatData *p = (IsoMatData *)tab;

	// ���e���W���E�Ð����E�΍����͍s��E�������͂̌v�Z
	G = p->young / (1 + p->poisson);
	calc_deviatoric_stres (estres, &sb, stres);
	sigbar = mieses_stres (stres);

	f=p->young/(1-p->poisson*p->poisson)/2;
	S[1]=f*(estres[0]+p->poisson*estres[1]);
	S[2]=f*(p->poisson*estres[0]+estres[1]);
	S[4]=2*G*estres[3];
	S[0]=4*sigbar*sigbar*hard/9+S[1]*estres[0]+S[2]*estres[1]+2*S[4]*estres[4];

	D[0][0]=-S[1]*S[1]/S[0];
	D[1][1]=-S[2]*S[2]/S[0];
	D[0][1]=D[1][0]=-S[1]*S[2]/S[0];
	D[0][2]=D[2][0]=-S[1]*S[4]/S[0];
	D[1][2]=D[2][1]=-S[2]*S[4]/S[0];
	D[2][2]=-S[4]*S[4]/S[0];
}

void calc_PSN_IsoDp_Mat (double *_D, const MatData *tab, const double *stres, double hard)
{
	double *D[3] = { _D, _D+3, _D+6, };
	double S[7], G, f, sb, sigbar, estres[6];
	IsoMatData *p = (IsoMatData *)tab;

	// ���e���W���E�Ð����E�΍����͍s��E�������͂̌v�Z
	G = p->young / (1 + p->poisson);
	calc_deviatoric_stres (estres, &sb, stres);
	sigbar = mieses_stres (stres);

	f = p->young / (1 - 2 * p->poisson) / (1 + p->poisson);
	S[1]=f*((1-p->poisson)*estres[1]+p->poisson*(estres[2]+estres[3]));
	S[2]=f*((1-p->poisson)*estres[2]+p->poisson*(estres[3]+estres[1]));
	S[3]=f*((1-p->poisson)*estres[3]+p->poisson*(estres[1]+estres[2]));
	S[4]=2*G*estres[4];
	S[0]=4*sigbar*sigbar*hard/9+S[1]*estres[1]+S[2]*estres[2]+S[3]*estres[3]+2*S[4]*estres[4];

	D[0][0]=-S[1]*S[1]/S[0];
	D[1][1]=-S[2]*S[2]/S[0];
	D[0][1]=-S[1]*S[2]/S[0]; D[1][0]=D[0][1];
	D[0][2]=-S[1]*S[4]/S[0]; D[2][0]=D[0][2];
	D[1][2]=-S[2]*S[4]/S[0]; D[2][1]=D[1][2];
	D[2][2]=-S[4]*S[4]/S[0];
}

void calc_AXSol_IsoDp_Mat (double *_D, const MatData *tab, const double *stres, double hard)
{
	int i, j;
	double *D[4] = { _D, _D+4, _D+8, _D+12, };
	double S, G, sb, sigbar, estres[6];
	IsoMatData *p = (IsoMatData *)tab;

	// ���e���W���E�Ð����E�΍����͍s��E�������͂̌v�Z
	G = p->young / (1 + p->poisson);
	calc_deviatoric_stres (estres, &sb, stres);
	sigbar = mieses_stres (stres);
	S = 2 * sigbar * sigbar / 3 * (1 + hard / G / 3);

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			D[i][j] = -estres[i] * estres[j] / S;
}

void calc_3DSolid_IsoDp_matrix (double *_D, const MatData *tab, const double *stres, double hard)
{
	int i, j;
	double *D[6] = { _D, _D+6, _D+12, _D+18, _D+24, _D+30, };
	double S, G, sb, sigbar, estres[6];
	IsoMatData *p = (IsoMatData *)tab;

	// ���e���W���E�Ð����E�΍����͍s��E�������͂̌v�Z
	G = p->young / (1 + p->poisson);
	calc_deviatoric_stres (estres, &sb, stres);
	sigbar = mieses_stres (stres);
	S = 2 * sigbar * sigbar / 3 * (1 + hard / G / 3);

	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++)
			D[i][j] = -estres[i] * estres[j] / S;
}

/*
 * �w���͂̒�`
 * �Ƃ肠�����ÓT�I��Plager��Ziegler�̒�`�ɂ����̂��f�ڂ���
 * ���l�^�͓n�Ӎ_�񎁂�FEM�{(�l�b�g�ł�PDF)�ɂ��B
 * ���ɂ���E��̃_�b�V���͕΍������ł���
 *
 * �p�����[�^�[
 *   alp_out: ���݂̔w���͑��x�B���[�`�����Ōv�Z
 *   stress:  ��������
 *   alp_in:  �����w����
 *   sig_par: ���͑��x
 */
// Plager�ɂ��w���͂̒�`�B
// ��ij,t=(��'ij-��'ij)(��'kl-��'kl)��kl,t/(2��y^2)
void calcPlagerBackstresFactor(double *alp_out, const double *stres, const double *alp_in, const double *sig_par)
{
	int i;
	double ea [6], es [6], sigy, sb, eb;
	calc_deviatoric_stres (es, &sb, stres);
	calc_deviatoric_stres (ea, &eb, alp_in);
	sigy=mieses_stres(stres);

	for (i = 0; i < 6; i++) {
		double f = 0;
		int j;
		for (j = 0; j < 3; j++)
			f += (es [j] - ea [j]) * sig_par [j];
		for (; j < 6; j++)
			f += 2 * (es [j] - ea [j]) * sig_par [j];
		alp_out [i] = (es [i] - ea [i]) * f / (2 * sigy * sigy);
	}
}

// Ziegler�ɂ��w���͂̒�`
// ��ij,t=(��ij-��ij)(��'kl-��'kl)��kl,t/(2��y^2)
void calcZieglerBackstresFactor(double *alp_out, const double *stres, const double *alp_in, const double *sig_par)
{
	int i;
	double ea [6], es [6], sigy, sb, eb;
	calc_deviatoric_stres (es, &sb, stres);
	calc_deviatoric_stres (ea, &eb, alp_in);
	sigy = mieses_stres (stres);

	for (i = 0; i < 6; i++) {
		double f = 0;
		int j;
		for (j = 0; j < 3; j++)
			f += (es [j] - ea [j]) * sig_par [j];
		for (; j < 6; j++)
			f += 2 * (es [j] - ea [j]) * sig_par [j];
		alp_out [i] = (stres [i] - alp_in [i]) * f / (2 * sigy * sigy);
	}
}

/* ���R�r�A���̋t�s����v�Z����B�Ƃ��ɍs��l��Ԃ��B */
double invJ_2D (double *iJ, const double *J)
{
	double detJ;
	detJ = J [0] * J [3] - J [2] * J [1];
	iJ [0] =  J [3] / detJ;
	iJ [1] = -J [1] / detJ;
	iJ [2] = -J [2] / detJ;
	iJ [3] =  J [0] / detJ;
	return detJ;
}

double invJ_3D (double *iJ, const double *J)
{
	double detJ;
	detJ =    J [0] * (J [4] * J [8] - J [5] * J [7])
		+ J [1] * (J [5] * J [6] - J [3] * J [8])
		+ J [2] * (J [3] * J [7] - J [4] * J [6]);
	iJ [0] = (J [4] * J [8] - J [5] * J [7]) /detJ;
	iJ [1] = (J [7] * J [2] - J [8] * J [1]) /detJ;
	iJ [2] = (J [1] * J [5] - J [2] * J [4]) /detJ;
	iJ [3] = (J [5] * J [6] - J [3] * J [8]) /detJ;
	iJ [4] = (J [8] * J [0] - J [6] * J [2]) /detJ;
	iJ [5] = (J [2] * J [3] - J [0] * J [5]) /detJ;
	iJ [6] = (J [3] * J [7] - J [4] * J [6]) /detJ;
	iJ [7] = (J [6] * J [1] - J [7] * J [0]) /detJ;
	iJ [8] = (J [0] * J [4] - J [1] * J [3]) /detJ;
	return detJ;
}

// MODULE B_MAT
/*
 * �v�f�����̌v�Z
 * �悭���鏑�ЂŐ�������Ă�����̂ƈقȂ�ȉ��̏����Ōv�Z���s��
 * �z�Ȍ`�ŗv�f���������߂Ă��Ȃ��B
 * ���R�r�A���̌v�Z
 *   �`��֐��̔����}�g���b�N�XH�Ɨv�f�̐ߓ_���W�l���������邱�Ƃŋ��߂���B
 *   ���R�����A�������o���V�Čv�Z����ƕێ�R�X�g�̑�����������Ƃ�����̂Œ��ӁB
 * ���R�r�A���̍s��l�Ƌt�s������߂�B
 * J^-1��H�̐ς����߂�B�e����͂Ȃ�A����łЂ��݁\�ψʍs������߂����ƂɂȂ�B
 * J^-1�EH�̊e��ɂ��čs��ς����߂�B
 * ���̌��ʂɕ����l�Ɛ��l�ϕ��̏d�݂����������̂����߂�B
 * �S�̍����s��ɑg�ݍ��ށB
 */

/* ���R�r�A���̌v�Z �`��֐��̔����`H�ɗv�f�̐ߓ_���W�l�������邱�Ƃœ�����B */
void calc_Jmat (double *J, const double *H, const Element *e)
{
	int i, j, k;
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	memset (J, 0, sizeof (double) * ndim * ndim);
	for (i = 0; i < ndim; i++)
		for (j = 0; j < ndim; j++)
			for (k = 0; k < nnode; k++)
				J [i * ndim + j] += H [i * nnode + k] * node [e->conn [k]].coord [j];
}

/*
   �񎟌��z��̊|���Z
   A[sz1][sz2] = B[sz1][sz3] * C [sz3][sz2]
 */
void mat_mul2 (double *A, const double *B, const double *C, int sz1, int sz2, int sz3)
{
	int i, j, k;
	memset (A, 0, sizeof (double) * sz1 * sz3);
	for (i = 0; i< sz1; i++)
		for (j = 0; j < sz2; j++)
			for (k = 0; k < sz3; k++)
				A [i * sz2 + j] += B [i * sz3 + k] * C [k * sz2 + j];
}

/*
   �v�f�����̌v�Z���@
   ����[1980]�ɂ��ƕψʂЂ��݃}�g���b�N�X��2�������ʉ���/�Ђ��ݗv�f�̏ꍇ�A�ȉ��̂悤�ɂȂ�B
   B = [ 1 0 0 0 ][ J^-1�EH ��      ]
       [ 0 0 0 1 ][ ��      J^-1�EH ]
       [ 0 1 1 0 ]
   ������J^-1�EH�ƃ���nnode�~ndim�̍s��ł���A���̓[���s��AJ^-1�̓��R�r�A���̋t�s��A
   H�͌`��֐��̔����`���}�g���b�N�X�z�u�������̂Ŋe��Ɍ`��֐���}��ϐ��ŕΔ����������̂����ԁB
   ������J^-1�EH���ȉ��̂悤�ɕ\������B
   J^-1�EH = [ H1 ]
             [ H2 ]
   �����B��W�J����ƈȉ��̂悤�ɂȂ�B

   B = [ 1 0 0 0 ][ H1 ��  ]
       [ 0 0 0 1 ][ H2 ��  ]
       [ 0 1 1 0 ][ �� H1 ]
                  [ �� H2 ]
     = [ H1 �� ]
       [ �� H2 ]
       [ H2 H1 ]
*/
void calc_PL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
{
	int i;
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	ElementCalcInfo *cp = e->cinfo + nint;
	int size = ndim * nnode;
	double J [ndim * ndim], iJ [ndim * ndim];
	calc_Jmat (J, H, e);
	cp->detJ = invJ_2D (iJ, J);
	mat_mul2 (cp->iJH, iJ, H, ndim, nnode, ndim);
	// B�}�g���b�N�X�̍쐬�B���΂�(�܂莩�R�x��)�ɒl�������Ă���B
	for (i = 0; i < nnode; i++) {
		B [0 * size + i * 2 + 0] = B [2 * size + i * 2 + 1] = cp->iJH [i];
		B [2 * size + i * 2 + 0] = B [1 * size + i * 2 + 1] = cp->iJH [i + nnode];
		B [1 * size + i * 2 + 0] = B [0 * size + i * 2 + 1] = 0;
	}
}

/*
   ���Ώ̌ő̗v�f�̍������v�Z����B
   ��{�I�ɓ񎟌����ʗv�f�Ɠ������W�b�N�����A�p�x�����̂Ђ��݂����邽�߂��̕����̌v�Z�ɑ��Ⴊ������B
   �Ђ��ݕψʍs��͎��̂悤�ɂȂ�B�a�����̂Ђ��݂�
       1/r=N/x
   �œ�����B������x�͐ϕ��_�ɂ�����e�ϕ��_�ł̑��a�ɂȂ�B
   B = [ 1 0 0 0 0 ][ J^-1*H 0       ]
       [ 0 0 0 1 0 ][ 0      J^-1*H  ]
       [ 0 0 0 0 1 ][ 1/r    0       ]
       [ 0 1 1 0 0 ]
   ������H3=1/x�ƒu���A���ʗv�f���l�ɓW�J����ƈȉ��̂悤�ɂȂ�B
   B = [ 1 0 0 0 0 ][ H1 0  ]
       [ 0 0 0 1 0 ][ H2 0  ]
       [ 0 0 0 0 1 ][ 0  H1 ]
       [ 0 1 1 0 0 ][ 0  H2 ]
                    [ H3 0  ]
     = [ H1 0  ]
       [ 0  H2 ]
       [ H3 0  ]
       [ H2 H1 ]
*/
void calc_AXSOL_B_mat (double *B, const Element *e, const double *N, const double *H, int nint)
{
	int i;
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	int size = ndim * nnode;
	ElementCalcInfo *cp = e->cinfo + nint;
	double J [ndim * ndim], iJ [ndim * ndim];
	double cx;
	calc_Jmat (J, H, e);
	cp->detJ = invJ_2D (iJ, J);
	mat_mul2 (cp->iJH, iJ, H, ndim, nnode, ndim);
	// �ϕ��_�ł̔��a�����߂�B
	cx = 0;
	for (i = 0; i < nnode; i++)
		cx += N [i] * node [e->conn [i]].coord [0];
	cp->iJH [nnode * 3] = cx;
	for (i = 0; i < nnode; i++) {
		B [0 * size + i * 2 + 0] = B [3 * size + i * 2 + 1] = cp->iJH [i];
		B [3 * size + i * 2 + 0] = B [1 * size + i * 2 + 1] = cp->iJH [i + nnode];
		B [2 * size + i * 2 + 0]                            = cp->iJH [i + nnode * 2] = N [i] / cx;
		B [1 * size + i * 2 + 0] = B [0 * size + i * 2 + 1] = B [2 * size + i * 2 + 1] = 0;
	}
}

/*
   �S�̍����̃Z�b�g�A�b�v
 */
// MODULE TOTAL_MAT_MISC

// �Ђ��ݕψʍs��Ɨv�f�����̐錾�B���\�ł����̂ł����Ő錾����B
double __DB__ [6 * 64 * sizeof (double)];
double __B__ [6 * 64 * sizeof (double)];

/*
 * �v�f�����}�g���b�N�X��S�̍����ɌJ�荞��
 * �e�ϕ��_���Ƃɍs���B
 */
void attach_B_mat (const Element *e, const double *D, int nint)
{
	int i, j, k, kx, ky, ypos, ex, ey, dofx, dofy;
	double Ke, thick;
	ElemIntegralInfo *info = e->info->info1;
	int nvstat = etinfo [prob].nvstat;
	int nnode = info->nnode;
	int ndof = etinfo [prob].ndof;
	int size = nnode * ndof;
	ElementCalcInfo *cp = e->cinfo + nint;
	// �v�f�����ɕK�v�Ȓl�̌v�Z�B�ϕ��_�ł̏d�݂ƃ��R�r�A���̐ς����߂�B
	double factor = info->wt [nint] * cp->detJ;
	// B^T*D���v�Z����B
	memset (__DB__, 0, size * nvstat * sizeof (double));
	mat_mul2 (__DB__, D, __B__, nvstat, size, nvstat);
	// �v�f�����̌v�Z�B�e�ʂ��傫�����߈�x�ɑS�����v�Z�������v�Z����B
	for (i = 0; i < size; i++) {
		ex = e->conn [i / ndof]; dofx = i % ndof;
		kx = rev_line_info [ex * mdof + dofx];
		if (kx >= rank_line_info [KLDOF_NONZERO][0]) continue;
		for (j = 0; j < size; j++) {
			ey = e->conn [j / ndof]; dofy = j % ndof;
			ypos = ey * mdof + dofy;
			ky = rev_line_info [ypos];
			// ������Ke[i][j]�̂݌v�Z����
			Ke = 0;
			for (k = 0; k < nvstat; k++)
				Ke += __B__ [k * size + i] * __DB__ [k * size + j];
			// �ŏI�I�ȗv�f�����̌v�Z�B���ʉ��́E�c���Ǝ��Ώ̖��ł͂�����W�����قȂ�B
			// �W���͐ϕ��_���Ƃ̏d�݂ƃ��R�r�A���̐ς����A���Ώ̖��ł͔��a��������B
			switch (prob) {
			case PT_PSS:
				thick = ((PSS_GeomData *) (e->part->geom))->thick;
				Ke *= thick;
			case PT_PSN:
				Ke *= factor;
				break;
			case PT_AXSOL:
				Ke *= cp->iJH [nnode * 3] * factor;
				break;
			}
			// �v�f������S�̍����ɑg�ݍ���
			// �o���h�}�g���b�N�X�v�Z�܂ł́u���̂܂܁v�ʒu��T���đg�ݍ��ނ��A
			// �Ώ̃}�g���b�N�X�͑Ίp���̏ゾ���ɒl�����邽�߁A�ʒu�𒲂ׂĊY������l�������Z�b�g����B
			if (ky < rank_line_info [KLDOF_NONE][1]) {
				if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END) {
					line_info [kx].p [ky] += Ke;
				} else if (sys_range_mode == RANGE_SYMM_UPPER) {
					if (ky >= kx)
						line_info [kx].p [ky] += Ke;
				}
			// �ψʍS�������̂����鎩�R�x�́A�v�f�����ƕψʂ����������̂��׏d�x�N�g���֑������ށB
			} else if (ky < rank_line_info [KLDOF_NONZERO][1])
				rhs_value [kx] -= Ke * disp [ypos];
		}
	}
}

/* ���E�����̍쐬
   �׏d�͌v�Z�p�z��ɓ����
 */
/* �׏d�f�[�^�[�̌v�Z */
static void calc_load_value (const LoadData *p, double *load)
{
	double v;
	int i, ndim;
	ndim = etinfo [prob].ndim;
	v = 1;
	switch (p->val_type % 10) {
	case BC_VAL:
	case BC_NORMAL:
		for (i = 0; i < ndim; i++)
			load [i] = p->value [i];
		break;
	case BC_VEC:
		for (i = 0, v = 0; i < ndim; i++)
			v += p->value [i] * p->value [i];
		v = sqrt (v);
		// 0�x�N�g���Ȃ�l�̐ݒ���������[�`�����甲����
		if (fabs (v) <= tolerance) break;
	case BC_VAL_VEC: // fall thrurh
		for (i = 0; i < ndim; i++)
			load [i] = p->value [ndim] * p->value [i] / v;
		break;
	}
}

// �ߓ_�׏d�̐ݒ�
void set_pload (const LoadData *p)
{
	int i, j, pos, rpos;
	double f;
	for (i = 0; i < mdof; i++) {
		if (fabs (p->load [i]) < tolerance) continue;
		for (j = 0; j < p->ndata; j++) {
			pos = p->data [j] * mdof + i;
			rpos = rev_line_info [pos];
			if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
			f = p->load [i];
			// ���Ώ̖��̏ꍇ�A�׏d�͕K�����a���|���邽�߁A�׏d�_��X���W�����o���悤�ɂ��Ă���
			if (prob == PT_AXSOL) f *= coord [p->data [j] * mdof];
			lhs_value [rpos] = rhs_value [rpos] += f;
		}
	}
}

/* �Ӊ׏d�̐ݒ�*/
static void calc_edge_load (LoadData *p, const Element *e, int *face, int nint)
{
	int i, j, k;
	ElemIntegralInfo *info = e->info->info3;
	int ndim = etinfo [prob].ndim, ndof = etinfo [prob].ndof;
	int nnode = info->nnode;
	int pos, rpos;
	double f, J, x [ndim], n [ndim];
	double *H = info->dN + nint * nnode;
	double *N = info->N + nint * nnode;
	double wt = info->wt [nint];
	// �ڐ��x�N�g���ƃ��R�r�A�������߂�
	memset (x, 0, sizeof (double) * ndim);
	for (J = i = 0; i < ndim; i++) {
		for (j = 0; j < nnode; j++) {
			x [i] += H [j] * node [e->conn [face [j]]].coord [i];
		}
		J += x [i] * x [i];
	}
	J = sqrt (J);
	if (fabs(J) < tolerance) return;
	for (i = 0; i < ndim; i++) x [i] /= J;	// �P�ʃx�N�g���ɂ���
	if (p->val_type % 10 == BC_NORMAL) {		// �ʂɐ����ȉ׏d�̌��������߂�
		if (ndim == 2)
			n [0] = x [1], n [1] = -x [0];
		if (ndim == 3)
			;	// ������
		for (i = 0; i < nnode; i++) {
			for (j = 0; j < ndof; j++) {
				pos = e->conn [face [i]] * ndof + j;
				rpos = rev_line_info [pos];
				if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
				f = J * wt * N [i] * (p->load [0] * n [j] + p->load [1] * x [j]);
				if (prob == PT_AXSOL) f *= x [0] * J;
				lhs_value [rpos] = rhs_value [rpos] += f;
			}
		}
	} else {
		for (i = 0; i < nnode; i++) {
			for (j = 0; j < ndof; j++) {
				pos = e->conn [face [i]] * ndof + j;
				rpos = rev_line_info [pos];
				if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
				f = J * wt * N [i] * p->load [j];
				if (prob == PT_AXSOL) f *= x [0] * J;
				lhs_value [rpos] = rhs_value [rpos] += f;
			}
		}
	}
}

static void set_fload (LoadData *p)
{
	int i, j, nnode, *elist;
	Element *e;
	for (i = 0; i < p->ndata; i++) {
		e = elem + p->data [i];
		nnode = e->info->info3->nnode;
		elist = e->info->edge_list + nnode * p->data [i * 2 + 1];
		for (j = 0; j < e->info->info3->ntint; j++)
			calc_edge_load (p, e, elist, j);
	}
}

// �ϕ��_���Ƃ̑̐ϗ͂����߂�
static void calc_bload (const LoadData *p, const Element *e, int nint)
{
	int i, j;
	ElemIntegralInfo *info = e->info->info4;
	int ndim = etinfo [prob].ndim, nnode = info->nnode;
	int ndof = etinfo [prob].ndof;
	int pos, rpos;
	double f, detJ, cx, J [ndim * ndim], iJ [ndim * ndim];
	double *H = info->dN + nint * nnode * ndim;
	double *N = info->N + nint * nnode;
	double wt = info->wt [nint];
	// ���R�r�A�������߂�B
	calc_Jmat (J, H, e);
	detJ = invJ_2D (iJ, J);
	f = detJ * wt;
	// ���d���l���B
	if (p->val_type / 10 == BC_DENSITY) f *= e->part->mat->density;
	// ���Ώ̖��ł͐ϕ��_�̔��a��������K�v������B
	if (prob == PT_AXSOL) {
		cx = 0;
		for (j = 0; j < nnode; j++)	// �ϕ��_�ł̔��a�����߂�
			cx += N [j] * coord [e->conn [j] * mdof];
		f *= cx;
	}
	for (j = 0; j < nnode; j++)
		for (i = 0; i < ndof; i++) {
			pos = e->conn [j] * ndof + i;
			rpos = rev_line_info [pos];
			if (rpos >= rank_line_info [KLDOF_NONE][1]) continue;
			if (fabs (p->load [i]) < tolerance) continue;
			lhs_value [rpos] = rhs_value [rpos] += p->load [i] * N [j] * f;
		}
}

static void set_bload (const LoadData *p)
{
	int i, j;
	for (i = 0; i < p->ndata; i++) {
		Element *e = &elem [p->data [i]];
		for (j = 0; j < e->info->info4->nnode; j++)
			calc_bload (p, e, j);
	}
}

static void calc_load ()
{
	int i;
	LoadData *p;
	for (i = 0, p = pload; i < npload; i++, p++) {
		calc_load_value (p, p->load);
		set_pload (p);
	}
	for (i = 0, p = fload; i < nfload; i++, p++) {
		calc_load_value (p, p->load);
		set_fload (p);
	}
	for (i = 0, p = bload; i < nbload; i++, p++) {
		calc_load_value (p, p->load);
		set_bload (p);
	}
}

// MODULE TOTAL_MAT
/*
   ��������S�̍����̍쐬�֘A���[�`���{��

   �ŏI�I�ɋ��߂�S�̍����s��͈ȉ��̂Ƃ���ł���B
   	u1:���m�ψ�
	u2:0�ł͂Ȃ����m�ψ�
	u3:0�̊��m�ψ�
	u4:�v�f�ɐڑ����Ă��Ȃ��ψ�
   [ K11 K12 K13 K14 ][ u1 ] = [ f1 ]
   [ K21 K22 K23 K24 ][ u2 ]   [ f2 ]
   [ K31 K32 K33 K34 ][ u3 ]   [ f3 ]
   [ K41 K42 K43 K44 ][ u4 ]   [ f4 ]

   ������W�J����ƈȉ��̂悤�ɂȂ�B
   K11*u1+K12*u2+K13*u3+K14*u4 = f1
   K21*u1+K22*u2+K23*u3+K24*u4 = f2
   K31*u1+K32*u2+K33*u3+K34*u4 = f3
   K41*u1+K42*u2+K43*u3+K44*u4 = f4

   u4�̍��͑Ή�����v�f�������������������Ȃ�����0�ƂȂ�B�]���ĕ���������֘A���鍀�����O�ł���B
   u3=0������A���̍����ق��ď��O�ł���B�������A�v�Z�㍄����0�ƂȂ�ƌv�Z�Ɏ��s����̂�
   ���͂���������̂ŁA���̓_���l�������K33=-I�Ƃ���K�v������B
   ����u2���E�ӂɈړ�����BK22��K33���lK22=-I�ƒu�����Ƃ��ł���B
   �W�J���͈ȉ��̂悤�ɒu����������B
   K11*u1     = f1-K12*u2
   K21*u1-If2 =   -K22*u2
   K31*u1-If3 =   -K32*u2

   ������܂Ƃ߂�ƈȉ��̂悤�ɂȂ�Bf2��f3�͖��m��(����)�ɂȂ�̂ŁA���ӂɏW�߂Ă���B
   [ K11 O   O   ][ u1 ] = [ f1-K12*u2 ]
   [ K21 -I  O   ][ f2 ]   [ -K22*u2 ]
   [ K31 O   -I  ][ f3 ]   [ -K32*u2 ]
   ��������̂܂܉����Ă������̂����A�Ђ��݂Ɖ��͂����߂���A�v�f�̓����ߓ_�͂����߂邱�Ƃ��ł��A
   ���̑��a�͔��͂ƊO�̘͂a�ɓ������Ȃ邱�Ƃ��m���Ă���B�ł���Ύ����I��1�s�ڂɂ��ĉ����΂������ƂɂȂ�B

   ���ɔ����v�Z���s���ꍇ���l����B
   ���R�����A���������v�f�������Čv�Z����͓̂���ł͂Ȃ�����A�S�̍�������x�ۑ����邱�ƂɂȂ�B

   .1) �e�v�f�̊e�ߓ_�ɂ��čő�l�ƍŏ��l��T���B
   .2) 1)�̊����̒l���͈͊O�ɂ���΁A�����̒l���X�V����
   .3) 1)�`2)���e���R�x�ɐݒ肷��B
   .1) ���q�ߓ_��T���B��������X�g�̖����ɓ���Ă����B
   .3) rev_line_info�𑖍����Aline_info�̑����ʒu�Ɏ��R�x�Ɛߓ_�ԍ����L�^����B
   .1) ��S�����R�x�A0�ł͂Ȃ��S�����R�x�A0�S�����R�x�A���q�̏��ɕ��ׂ�B
   .2) line_info�͂��̏��Ԃɕ��ׂ�B�Ȃ��A�K�v�Ȃ͔̂�S�����R�x�����ł���B
*/

/* �v�f�ɐڑ����Ă��Ȃ��ߓ_��T�� */
void set_connected_dof ()
{
	int i, j, k, nid;
	int ndof = etinfo [prob].ndof;
	Element *p;

	// ���߂ɐ��q�ߓ_�ł���Ɖ��肷��B
	for (i = 0; i < ntnode * mdof; i++)
		nodal_dof_list [i] = KLDOF_ORPHAN;
	// �R�l�N�e�B�r�e�B�̃`�F�b�N�B�v�f�ɐڑ�����ߓ_��KLDOF_NONE���Z�b�g
	for (i = 0, p = elem; i < ntelem; i++, p++)
		for (j = 0; j < p->info->info1->nnode; j++) {
			nid = p->conn [j];
			for (k = 0; k < ndof; k++)
				nodal_dof_list [nid * mdof + k] = KLDOF_NONE;
		}
}

/* ���E������K�p���鎩�R�x�Ɉ������ */
void set_bc_pos ()
{
	int i, j, k, pos;
	FixDisp *fd;

	for (i = 0, fd = fdisp; i < nfdisp; i++, fd++) {
		for (j = 0; j < mdof; j++) {
			if (!fd->flags [j]) continue;
			for (k = 0; k < fd->ndata; k++) {
				int np = fd->data [k];
				pos = np * mdof + j;
				//
				if (nodal_dof_list [pos] != KLDOF_ORPHAN)
					nodal_dof_list [pos] = (fd->value [j] > 0)? KLDOF_NONZERO: KLDOF_ZERO;
				if (fd->value [j] > 0) disp [pos] += fd->value [j];
			}
		}
	}
}

/* ��S�����R�x�E�S�����R�x�E���q�ߓ_�̏��ɐߓ_�ԍ��Ǝ��R�x����ׂ� */
void create_rev_line_info ()
{
	int i, j, kltype, pos, count;
	// �ŏ��ƍŌ�̈ʒu�͖ق��Č��܂�̂ł��̎��_�Őݒ�
	rank_line_info [KLDOF_NONE][0] = 0;
	rank_line_info [KLDOF_ORPHAN][1] = ntnode * mdof;
	for (kltype = KLDOF_NONE, count = 0; kltype < KLDOF_ORPHAN; kltype++) {
		rank_line_info [kltype][1] = rank_line_info [kltype][0];
		// �S���^�C�v�̌����Ɣz�u�B�z���[1]-[0]�����ɂȂ�B
		// �]���Ċe�f�[�^�[�̏I�_�ʒu�͎��ۂ̏I�_�ʒu�̈�O�ɂ���B
		// ����ď����A�N�Z�X����ꍇ��C�̔z��ɕ키�悤�ɂ��Ă���B
		for (i = 0; i < ntnode; i++) {
			for (j = 0; j < mdof; j++) {
				pos = i * mdof + j;
				if (nodal_dof_list [pos] == kltype) {
					decode_rev_info [count] = pos;
					rev_line_info [pos] = count++;
				}
			}
		}
		// ���̃f�[�^�[�̋N�_�ʒu�͑O�̃f�[�^�[�́u�I�_�ʒu�v�ɓ������B
		rank_line_info [kltype + 1][0] = rank_line_info [kltype][1] = count;
	}
	rank_line_info [KLDOF_ORPHAN][0] = rank_line_info [KLDOF_ZERO][1];
}

void init_K_line_info ()
{
	int i, j;
	K_line_info *p;
	for (i = j = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		while (rev_line_info [j] >= rank_line_info [KLDOF_NONE][1])
			j++;
		p->start = p->end = i;
		j++;
	}
}

void set_range_all ()
{
	int i, size = rank_line_info [KLDOF_NONE][1] - rank_line_info [KLDOF_NONE][0];
	K_line_info *p;
	K_line_info tmp = {
		0,
		rank_line_info [KLDOF_NONE][0],
		rank_line_info [KLDOF_NONE][1],
		size, NULL,
	};

	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		*p = tmp;
		p->pos = size * i;
	}
}

/* �e�v�f�̃R�l�N�e�B�r�e�B�\����e���C���ɂ�����o���h�������߂�B */
void set_range_by_elem ()
{
	int i, j, k, pos, minpos, maxpos;
	int ndof = etinfo [prob].ndof;
	K_line_info *p, *pp;
	Element *e;

	for (i = 0, e = elem; i < ntelem; i++, e++) {
		ElemIntegralInfo *info = e->info->info1;
		int nnode = info->nnode;
		minpos = ntnode * mdof;
		maxpos = 0;
		// �v�f�ɂ�����R�l�N�e�B�r�e�B�͈̔͂�T��
		for (j = 0; j < nnode; j++) {
			for (k = 0; k < ndof; k++) {
				pos = rev_line_info [e->conn [j] * mdof + k];
				if (pos >= rank_line_info [KLDOF_NONE][1]) continue;
				if (minpos > pos) minpos = pos;
				if (maxpos < pos) maxpos = pos;
			}
		}
		for (j = 0; j < nnode; j++) {
			for (k = 0; k < ndof; k++) {
				pos = rev_line_info [e->conn [j] * mdof + k];
				if (pos >= rank_line_info [KLDOF_NONE][1]) continue;
				p = line_info + pos;
				if (sys_range_mode == RANGE_BEGIN_END && p->start > minpos)
					p->start = minpos;
				if (p->end < maxpos + 1) p->end = maxpos + 1;
			}
		}
	}
	line_info [0].size = line_info [0].end - line_info [0].start;
	for (i = 1, pp = line_info, p = pp + 1; i < rank_line_info [KLDOF_NONE][1]; i++, p++, pp++) {
		if (pp->end > p->end) p->end = pp->end;
		p->size = p->end - p->start;
	}
}

/* �S�̍����}�g���b�N�X�̃T�C�Y���v�Z���� */
int set_K_size ()
{
	int i, size;
	K_line_info *p;

	for (i = size = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		p->pos = size;
		size += p->size;
	}
	return size;
}

/*
 * �S�̍����̑g�ݗ���
 */

/*
 * �S�̍����̃Z�b�g�A�b�v
 *   �ȉ��Ɏg�p���������ŏ����ɂȂ�悤�ɂ��邩�l����
 * �Z�@
 *   �S���R�x�����߂�
 *   ���C�����̃������[�����蓖�Ă�
 *   ���C�����̐ݒ���s���B�����ʒu�͑Ίp����ɂ���A�S���͂Ȃ���ԂƂ���B
 */
void create_K_info (FILE *fout)
{
	int i, size;
	K_line_info *p;

	size = ntnode * mdof;
	nodal_dof_list = calloc (sizeof (int), size);
	disp = calloc (sizeof (double), size);
	force = calloc (sizeof (double), size);
	rev_line_info = calloc (sizeof (int), size);
	decode_rev_info = calloc (sizeof (int), size);

	set_connected_dof ();
	set_bc_pos ();
	create_rev_line_info ();
	dump_rank_info (OFF);
	dump_nodal_dof (OFF);
	dump_rev_line_info (OFF);

	line_info = calloc (sizeof (K_line_info), rank_line_info [KLDOF_NONE][1]);
	rhs_value = calloc (sizeof (double), rank_line_info [KLDOF_NONE][1]);
	lhs_value = calloc (sizeof (double), rank_line_info [KLDOF_NONE][1]);
	init_K_line_info ();

	if (sys_range_mode == RANGE_ALL)
		set_range_all ();
	else
		set_range_by_elem ();

	size = set_K_size ();
	fprintf (stderr, "Total matrix size: " LLU0 "bytes\n", size * sizeof (double));
	fprintf (fout, "\nTotal matrix size: " LLU0 "bytes\n", size * sizeof (double));

	sysk = calloc (sizeof (double), size);
	if (sysk == NULL) {
		fputs ("Error: Cannot allocate system stiff matricx\n", stderr);
		free_data ();
		exit (-11);
	}
	dump_pline (OFF, OFF);

	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++)
		p->p = sysk + p->pos - p->start;
	dump_pline (OFF, OFF);
	dump_pline (OFF, OFF);
}

typedef void (*FP_calc_B) (double *, const Element *, const double *, const double *, int);
void calc_K_mat (FILE *fout)
{
	int i, j;
	int nvstat = etinfo [prob].nvstat;
	FP_calc_B calc_B [] = { NULL, calc_PL_B_mat, calc_PL_B_mat, calc_AXSOL_B_mat, NULL, };
	create_K_info (fout);
	calc_load ();
	dump_bc (OFF);
	dump_rank_info (OFF);
	dump_nodal_dof (OFF);
	dump_rev_line_info (OFF);
	dump_bc (OFF);
	for (i = 0; i < ntelem; i++) {
		Element *e = elem + i;
		ElemIntegralInfo *info = e->info->info1;
		double *N = info->N;
		double *H = info->dN;
		double *D = e->part->init_mat;
		int ndim = etinfo [prob].ndim, nnode = info->nnode;
		for (j = 0; j < info->ntint; j++) {
			(*calc_B [prob]) (__B__, e, N, H, j);
			attach_B_mat (e, D, j);
			N += nnode;
			H += nnode * ndim;
			dump_mat2 (e->cinfo [j].iJH, ndim, nnode, OFF);
			dump_mat2 (__B__, nvstat, nnode * ndim, OFF);
		}
	}
	dump_pline (OFF, OFF);
}

/*
 * �o���h�}�g���b�N�X�`���̕�����������
 */
void solve_band_mat (double tol)
{
	int i, j, k, size;
	double pivot, Kji;

	// �O�i����
	size = rank_line_info [KLDOF_NONE][1];
	for (i = 0; i < size; i++) {
		pivot = line_info [i].p [i];
		if (fabs (pivot) < tol) {
			fprintf (stderr, "Error: pivot tolerance under %g at %d\n", pivot, i);
			free_data ();
			exit (-12);
		}
		for (j = i + 1; j < line_info [i].end; j++) {
			//if (line_info [j].start > i) continue;
			if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END)
				Kji = line_info [j].p [i];
			else if (sys_range_mode == RANGE_SYMM_UPPER)
				Kji = line_info [i].p [j];
			if (fabs (Kji) < tol) continue;
			// �t���}�g���b�N�X�ƃo���h�}�g���b�N�X�̓s�{�b�g�̒�������v�Z���邪
			// �㔼���͑Ίp�v�f����v�Z����悤�ɂ���
			if (sys_range_mode == RANGE_ALL || sys_range_mode == RANGE_BEGIN_END)
				for (k = i + 1; k < j; k++)
					line_info [j].p [k] -= line_info [i].p [k] * Kji / pivot;
			for (k = j; k < line_info [i].end; k++)
				line_info [j].p [k] -= line_info [i].p [k] * Kji / pivot;
			rhs_value [j] -= rhs_value [i] * Kji / pivot;
		}
	}
	//fputs ("end of forword deletion\n", stderr);
	// ��ޑ��
	for (i = size - 1; i >= 0; i--) {
		pivot = line_info [i].p [i];
		for (j = i + 1; j < line_info [i].end; j++) {
			Kji = line_info [i].p [j];
			if (fabs (Kji) < tol) continue;
			rhs_value [i] -= rhs_value [j] * Kji;
		}
		rhs_value [i] /= pivot;
		line_info [i].p [i] /= pivot;
		//fprintf (stderr, "%5d%15g%15g\n", i, rhs_value [i], pivot);
	}
	dump_bc (OFF);
}

/* ���͂ƂЂ��݁A�ό`���z�e���\�����v�Z���� */
/*
 * �����Ђ��݂̌v�Z
 * �����Ђ��݂̌v�Z�̓�=Bu�Őϕ��_���Ƃɋ��܂�B
 * B�͕��ʉ��́E���ʂЂ��ݗv�f�̏ꍇ�ȉ��̂悤�ɂȂ�
   B = [ 1 0 0 0 ][ J^-1�EH ��      ]
       [ 0 0 0 1 ][ ��      J^-1�EH ]
       [ 0 1 1 0 ]
   �����Ń���nnode�~ndim�̃[���s��ł���BJ^-1�̓��R�r�A���̋t�s��A
   H�͌`��֐��̔����`���}�g���b�N�X�z�u�������̂Ŋe��Ɍ`��֐���}��ϐ��ŕΔ����������̂����ԁB
   ������J^-1�EH���ȉ��̂悤�ɕ\������B
   J^-1�EH = [ H1 ]
             [ H2 ]
   �����B��W�J����ƈȉ��̂悤�ɂȂ�B

   B = [ 1 0 0 0 ][ H1 �� ]
       [ 0 0 0 1 ][ H2 �� ]
       [ 0 1 1 0 ][ �� H1 ]
                  [ �� H2 ]
     = [ H1 �� ]
       [ �� H2 ]
       [ H2 H1 ]
   �c�Â̓�=B*u�ŋ��߂���B
   �� = [ H1 �� ][ u ]
        [ �� H2 ][ v ]
        [ H2 H1 ]
      = [ H1*u      ]
        [ H2*v      ]
        [ H2*u+H1*v ]

 * ���͂̌v�Z
 * ���͂ƂЂ��݂̊֌W�� ��=D�� ����e�ՂɌv�Z�����B
 *
 * �����ߓ_�͂̌v�Z
 * ���̓Ђ�p���� pe=��B^T��|J|dt1dt2 ����v�Z����B
 *
 * �ό`���z�̌v�Z
 *   �ȉ��̎�����v�Z���s���B
 *   �@F=(iJH*x)^T
 *   ������iJH�̓��R�r�A���̋t�s��ƌ`��֐��̌��z�Ax�͕ό`��̍��W�ł���B
 *   ���̂܂܋��߂�Ɠ]�u�`���ɂȂ�̂ŁA�v�Z����F�̍s������ւ���
 *
 * ���Ώ̌ő̗v�f�̌��ʂ��v�Z����B
   ��{�I�ɓ񎟌����ʗv�f�Ɠ������W�b�N�����A�a�����̂Ђ��݂����邽�߂��̕����̌v�Z�ɑ��Ⴊ������B
   �Ђ��ݕψʍs��͎��̂悤�ɂȂ�B
   B = [ 1 0 0 0 0 ][ J^-1*H 0       ]
       [ 0 0 0 1 0 ][ 0      J^-1*H  ]
       [ 0 0 0 0 1 ][ 0      dNi/dNj ]
       [ 0 1 1 0 0 ]
   �����Ōa�����̂Ђ��݂�dNi/dNj�ƒ�`�������A����͐ϕ��_�̒l�̑��a���Ƃ�A
   ���q�͒��ڂ���ϕ��_�ł̒l���Ƃ�B
   ������dR=dNi/dNj�ƒu���A���ʗv�f���l�ɓW�J����ƈȉ��̂悤�ɂȂ�B
   B = [ 1 0 0 0 0 ][ H1 0  ]
       [ 0 0 0 1 0 ][ H2 0  ]
       [ 0 0 0 0 1 ][ 0  H1 ]
       [ 0 1 1 0 0 ][ 0  H2 ]
                    [ 0  dR ]
     = [ H1 0  ]
       [ 0  H2 ]
       [ 0  dR ]
       [ H2 H1 ]
�@ ����ɕψʂ�������ƂЂ��݂����܂�B
   ��= [ ��x ��y ��z ��xy ]^T
   Bu = [ H1 0  ][ u ]
        [ 0  H2 ][ v ]
        [ 0  dR ]
        [ H2 H1 ]
      = [ H1*u      ]
        [ H2*v      ]
	[ dR*v      ]
	[ H2*u+H1*v ]
 */
void calc_PL_stress_strain ()
{
	int i, j, k, l, m, nnode, ntint, ndof, nvstat, ndim;
	double f, *D, *wt, *H1, *H2, *H3, *stress, *strain;
	Element *e;
	ElementCalcInfo *cp;

	nvstat = etinfo [prob].nvstat;
	ndof = etinfo [prob].ndof;
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		ElemIntegralInfo *info = e->info->info1;
		nnode = info->nnode; ntint = info->ntint;
		ndim = etinfo [prob].ndim;
		wt = info->wt;
		IsoMatData *mat = (IsoMatData *)(e->part->mat);
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			H1 = cp->iJH; H2 = H1 + nnode; H3 = H2 + nnode;
			strain = cp->strain; stress = cp->stress;
			memset (strain, 0, sizeof (double) * 6);
			memset (stress, 0, sizeof (double) * 6);
			//fprintf (stderr, "%5d\n", prob);
			// ���͂ƂЂ��݂̌v�Z�B
			for (k = 0; k < nnode; k++) {
				int pos1 = e->conn [k] * mdof;
				int pos2 = pos1 + 1;
				strain [0] += H1 [k] * disp [pos1];
				strain [1] += H2 [k] * disp [pos2];
				strain [3] += H2 [k] * disp [pos1] + H1 [k] * disp [pos2];
				if (prob == PT_AXSOL)
					strain [2] += H3 [k] * disp [pos1];
			}
			if (prob == PT_PSS || prob == PT_PSN)
				strain [2] = strain [3];
			D = e->part->init_mat;
			for (k = 0; k < nvstat; k++) {
				for (l = 0; l < nvstat; l++, D++)
					stress [k] += *D * strain [l];
			}
			if (prob == PT_PSS || prob == PT_PSN) {
				stress [3] = stress [2];
				switch (prob) {
					case PT_PSS:
						stress [2] = 0;
						strain [2] = mat->poisson / mat->young * (stress [0] + stress [1]);
						break;
					case PT_PSN:
						stress [2] = mat->poisson * (stress [0] + stress [1]);
						strain [2] = 0;
						break;
				}
			}
			// �ߓ_�͂̌v�Z�B�O�͂Ɣ��͂�������B
			f = cp->detJ * wt [j];
			for (k = 0; k < nnode; k++) {
				double f2 = f;
				int pos1 = e->conn [k] * ndof;
				int pos2 = pos1 + 1;
				if (prob == PT_AXSOL) f2 *= cp->iJH [nnode * 3];
				force [pos1] += f2 * (H1 [k] * stress [0] + H2 [k] * stress [3]);
				force [pos2] += f2 * (H2 [k] * stress [1] + H1 [k] * stress [3]);
				if (prob == PT_AXSOL) force [pos1] += f2 * (H3 [k] * stress [2]);
			}
			// �ό`���z�e���\���̌v�Z
			// F=(J^-1H)^T�ŋ��߂���B
			memset (cp->FT, 0, sizeof (double) * 9);
			for (k = 0; k < ndim; k++)
				for (l = 0; l < ndim; l++)
					for (m = 0; m < nnode; m++) {
						int node_id = e->conn [m];
						double xi = node [node_id].coord [l] + disp [node_id * ndim + l];
						cp->FT [l * 3 + k] += cp->iJH [k * nnode + m] * xi;
					}
			if (prob <= PT_AXSOL) cp->FT [8] = 1;
		}
	}
}

void calc_result ()
{
	int i, j, pos2, pos;
	// �ߓ_�ψʂƉ׏d�̐���
	for (i = 0; i < ntnode; i++) {
		for (j = 0; j < mdof; j++) {
			pos = i * mdof + j;
			if (pos < rank_line_info [KLDOF_NONE][1]) {
				pos2 = decode_rev_info [pos];
				disp [pos2] = rhs_value [pos];
			}
		}
	}
	calc_PL_stress_strain ();
}

/* ���ʏo�� */
void write_result (FILE *fout)
{
	int i, j, k, l;
	double *p;
	Element *e;
	ElementCalcInfo *cp;

	fprintf (fout, "Nodal displacement\n%10s%15s%15s%15s\n", "Number", "X-dir", "Y-dir", "Z-dir");
	for (i = 0, p = disp; i < ntnode; i++, p += mdof) {
		fprintf (fout, "%10d", node [i].number);
		for (j = 0; j < mdof; j++)
			fprintf (fout, "%15g", p [j]);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
	fprintf (fout, "Reaction/External force\n%10s%15s%15s%15s\n", "Number", "X-dir", "Y-dir", "Z-dir");
	for (i = 0, p = force; i < ntnode; i++, p += mdof) {
		fprintf (fout, "%10d", node [i].number);
		for (j = 0; j < mdof; j++)
			fprintf (fout, "%15g", p [j]);
		fputc ('\n', fout);
	}
	fputc ('\n', fout);
	fprintf (fout, "Small strain\n%10s%10s%15s%15s%15s%15s%15s%15s\n",
			"Number", "IPT", "X-dir", "Y-dir", "Z-dir", "XY-shar", "YZ-shar", "ZX-shar");
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			fprintf (fout, "%10d%10d", e->number, j + 1);
			for (k = 0; k < 6; k++)
				fprintf (fout, "%15g", cp->strain [k]);
			fputc ('\n', fout);
		}
	}
	fputc ('\n', fout);
	fprintf (fout, "Cauchy stress\n%10s%10s%15s%15s%15s%15s%15s%15s\n",
			"Number", "IPT", "X-dir", "Y-dir", "Z-dir", "XY-shar", "YZ-shar", "ZX-shar");
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			fprintf (fout, "%10d%10d", e->number, j + 1);
			for (k = 0; k < 6; k++)
				fprintf (fout, "%15g", cp->stress [k]);
			fputc ('\n', fout);
		}
	}
	fputc ('\n', fout);
	fprintf (fout, "Displacement gradient tensor\n%10s%10s%15d%15d%15d\n", "Number", "IPT", 1, 2, 3);
	for (i = 0, e = elem; i < ntelem; i++, e++) {
		int ntint = e->info->info1->ntint;
		for (j = 0, cp = e->cinfo; j < ntint; j++, cp++) {
			for (k = 0; k < 3; k++) {
				if (k == 0)
					fprintf (fout, "%10d%10d", e->number, j + 1);
				else
					fprintf (fout, "%20c", ' ');
				for (l = 0; l < 3; l++)
					fprintf (fout, "%15g", cp->FT [k * 3 + l]);
				fputc ('\n', fout);
			}
		}
	}
	fputc ('\n', fout);
}

/*
 * �������̉��
 */
/* index�܂ł̗v�f���f�[�^�[��S�폜���A�v�f���X�g���폜����B */
void free_elem_data (int index)
{
	int i, j;
	Element *p;
	ElementCalcInfo *cp;
	for (i = 0, p = elem; i <= index; i++, p++) {
		int ntint = p->info->info1->ntint;
		free (p->conn);
		for (j = 0, cp = p->cinfo; j < ntint; j++, cp++)
			free (cp->iJH);
		free (p->cinfo);
	}
}

void free_set_data (SetData *p)
{
	if (p == NULL) return;
	free (p->name);
	free (p->data);
}

void free_mesh_data ()
{
	int i;
	free (title); title = NULL;
	free (part); part = NULL;
	for (i = 0; i < nmat; i++) {
		if (mat [i] != NULL)
			free (mat [i]->name);
		free (mat [i]);
	}
	free (mat);
	nmat = 0; mat = NULL;
	for (i = 0; i < ngeom; i++) {
		if (geom [i] != NULL)
			free (geom [i]->name);
		free (geom [i]);
	}
	free (geom);
	ngeom = 0; geom = NULL;
	free (coord);
	free (node);
	ntnode = 0; node = NULL;
	free_elem_data (ntelem - 1);
	free (elem);
	ntelem = 0; elem = NULL;
	for (i = 0; i < nset; i++)
		free_set_data (set + i);
	free (set);
	nset = 0; set = NULL;
}

void free_load (int ndata, LoadData *p)
{
	int i;
	for (i = 0; i < ndata; i++)
		free_set_data ((SetData *)(p + i));
}

void free_bc_data ()
{
	int i;
	for (i = 0; i < nfdisp; i++)
		free_set_data ((SetData *)(fdisp + i));
	free (fdisp);
	nfdisp = 0; fdisp = NULL;
	free_load (npload, pload);
	free (pload);
	npload = 0; pload = NULL;
	free_load (nfload, fload);
	free (fload);
	nfload = 0; fload = NULL;
	free_load (nbload, bload);
	free (bload);
	nbload = 0; bload = NULL;
}

void free_sys_info ()
{
	free (line_info);
	free (nodal_dof_list);
	free (rev_line_info);
	free (decode_rev_info);
	free (sysk);
	free (rhs_value);
	free (lhs_value);
	free (disp);
	free (force);
}

void free_data ()
{
	free_mesh_data ();
	free_bc_data ();
	free_sys_info ();
}

int main (int argc, char **argv)
{
	FILE *fin, *fout;
	if (argc != 3) {
		fprintf (stderr, "usage: %s input output\n", argv [0]);
		exit(-1);
	}
	fprintf (stderr, "Start of analysis\n");
	init_einfo ();
	fprintf (stderr, "End of initializing\n");

	fin = fopen (argv [1], "r");
	if (fin == NULL) {
		fprintf (stderr, "cannot open input file: %s", argv [1]);
		return -1;
	}
	fout = fopen (argv [2], "w");
	if (fout == NULL) {
		fprintf (stderr, "cannot open output file: %s", argv [2]);
		return -1;
	}
	io_data (fin, fout);
	io_bc (fin, fout);
	fprintf (stderr, "End of reading data\n");

	calc_D_mat();
	sys_range_mode = RANGE_ALL;
	sys_range_mode = RANGE_BEGIN_END;
	sys_range_mode = RANGE_SYMM_UPPER;

	calc_K_mat (fout);
	fprintf (stderr, "End of construct stiffness matrix\n");
	solve_band_mat (tolerance);
	calc_result ();
	write_result (fout);
	fprintf (stderr, "End of solve probrem\n");

	free_data ();
	fprintf (stderr, "End of analysis\n");
	return 0;
}



