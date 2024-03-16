#include <c6x.h>
#include <math.h>
#include <stdlib.h>
#include "include\\var2.h"
#include "include\\myproject.h"
#include "include\\evmdm6437.h"
#include "include\\tracker.h"
#include "include\\IMG_sobel_3x3_8.h"
#include "include\\IMG_conv_5x5_i8_c8s.h"
#include "include\\IMG_conv_5x5_i8_c16s.h"
#include "include\\IMG_conv_3x3_i8_c8s.h"
#include "include\\ras_constant.h"


#pragma DATA_SECTION(m_morphology, "L1D_SRAM");
#pragma DATA_SECTION(d_morphology, "L1D_SRAM");

#pragma DATA_ALIGN(m_morphology, 32);
#pragma DATA_ALIGN(d_morphology, 32);

unsigned int m_morphology[10];
unsigned int d_morphology[10];

extern unsigned char IS1[2*512*36];
extern unsigned char ISD[6*512*36];
extern unsigned char ISD3[6*512*36];
extern unsigned char ISDResize[6*512*36];
extern unsigned char ISD3Resize[6*512*36];
double thd_b_old=0,thd_w_old=0;
void Binary_Image(unsigned char *restrict ID, int IHEIGHT, int IWIDTH,short hsn, short ws,short ys,short xs,unsigned char IR_Day_Select)
{

	int LengthIS,hs,i,xsn,ysn,i1,j,j2,LengthISn,hsd,hsn1,ws2,rem,l,ln,ld;
	int k,num_block;
	unsigned char *ISEdge,*ISEdge2,*ISEdge3;
	float	thd_w_f,thd_b_f,average_total;
	unsigned char	max_gray,min_gray;
	unsigned char	thd_w,thd_b;
	unsigned int th_b_3210,th_w_3210;
	int shift1=8;
	char mask[9]={1, 1, 1,
			1, 1, 1,
			1, 1, 1
			};
	char mask2[9]={0, 0, 0,
			1, 1, 1,
			0, 0, 0
			};
	short gaussian4_5x5[25]=
	//sigma=0.5   *256
	/*	{
	0 ,   0,     0 ,    0,     0,
	0 ,   3,    21,    3,    0,
	0 ,   21,    158,    21,    0,
	0 ,   3,    21,    3,    0,
	0 ,   0,     0 ,    0 ,    0*/
	 // sigma=0.75 *256
	{
		0 ,   1,     2 ,    1,     0, 
		1 ,   12,    30,    12,    1,
		2 ,   30,    72,    30,    2,
		1 ,   12,    30,    12,    1,
		0 ,   1,     2 ,    1 ,    0
	};

	m_morphology[0] = mask[0] > 0 ? ~0U : 0;
	m_morphology[1] = mask[1] > 0 ? ~0U : 0;
	m_morphology[2] = mask[2] > 0 ? ~0U : 0;
	m_morphology[3] = mask[3] > 0 ? ~0U : 0;
	m_morphology[4] = mask[4] > 0 ? ~0U : 0;
	m_morphology[5] = mask[5] > 0 ? ~0U : 0;
	m_morphology[6] = mask[6] > 0 ? ~0U : 0;
	m_morphology[7] = mask[7] > 0 ? ~0U : 0;
	m_morphology[8] = mask[8] > 0 ? ~0U : 0;
	d_morphology[0] = mask[0] <= 0 ? ~0U : 0;
	d_morphology[1] = mask[1] <= 0 ? ~0U : 0;
	d_morphology[2] = mask[2] <= 0 ? ~0U : 0;
	d_morphology[3] = mask[3] <= 0 ? ~0U : 0;
	d_morphology[4] = mask[4] <= 0 ? ~0U : 0;
	d_morphology[5] = mask[5] <= 0 ? ~0U : 0;
	d_morphology[6] = mask[6] <= 0 ? ~0U : 0;
	d_morphology[7] = mask[7] <= 0 ? ~0U : 0;
	d_morphology[8] = mask[8] <= 0 ? ~0U : 0;
	ysn=ys;
	xsn=xs;
	//ws*hsn must be multiple of 8, ws must be multiple of 64,hsn must be multiple of 6 and 2 ,
	hsn1=hsn;
	rem=hsn%48;
	if (rem == 0)
		hsn=hsn;
	else
		hsn=hsn+(48-hsn%48);
	num_block=6;
	hs=hsn/num_block;//+1;
	ws2=ws;
	rem=ws%64;
	if (rem == 0)
		ws=ws;
	else
		ws=ws+(64-ws%64);
	LengthISn=(hsn+6)*ws;
	LengthIS=(hs+1)*ws;
	ISEdge=(unsigned char*)malloc((3*LengthIS)*sizeof(char));
	ISEdge2=(unsigned char*)malloc((3*LengthIS)*sizeof(char));
	j2=0;
	i1=2*ws;
	l=ws*hs;
	ln=ws*hsn;
	for (j=0;j<2;j++)
	{
		DMA_fieldToRectOffset1DMA4ABSync( ID, IWIDTH*2, 2, ysn, xsn, hs, ws, IS1);
		DMA_Wait_total();
		RAS_UYVY_TO_Y(IS1, l, ISEdge+j2);
		j2=l+j2;
		ysn=ysn+hs;
		DMA_fieldToRectOffset1DMA4ABSync( ID, IWIDTH*2, 2, ysn, xsn, hs, ws, IS1);
		DMA_Wait_total();
		RAS_UYVY_TO_Y(IS1, l, ISEdge+j2);
		j2=l+j2;
		ysn=ysn+hs;
		DMA_fieldToRectOffset1DMA4ABSync( ID, IWIDTH*2, 2, ysn, xsn, hs, ws, IS1);
		DMA_Wait_total();
		RAS_UYVY_TO_Y(IS1, l, ISEdge+j2);
		j2=l+j2;
		ysn=ysn+hs;
		DMA_fieldToRectOffset1DMA4ABSync( ID, IWIDTH*2, 2, ysn, xsn, 2, ws, IS1);
		DMA_Wait_total();
		RAS_UYVY_TO_Y(IS1, 2*ws, ISEdge+j2);
		j2=0;
		ysn=ysn-2;
		Zero_Func_Mat(ISEdge2,3*LengthIS);
		for (i=0;i<(hsn/2+2);i++)
			IMG_conv_5x5_i8_c16s( ISEdge+i*ws, ISEdge2+i*ws+2,ws-4,ws, gaussian4_5x5,shift1);
		DMA_fieldToRectDMA4ABSync( ISEdge2, ws, 1,0,0, hsn/2, ws, ISD+i1);//hs must be multiple of 8
		DMA_Wait_total();
		DMA_fieldToRectDMA4ABSync( ISEdge, ws, 1,0,0, hsn/2, ws, ISD3+i1);//hs must be multiple of 8
		DMA_Wait_total();
		i1=ln/2;

	}
	free(ISEdge);
	free(ISEdge2);
    ISEdge=(unsigned char*)malloc(LengthISn*sizeof(char));
	DMA_fieldToRectDMA4ABSync( ISD, ws, 1,4,4,  hsn-8,ws-8, ISEdge);//ws must be multiple of 8
	DMA_Wait_total();
	int lnnew;
	lnnew=(hsn-8)*(ws-8);
	k=0;
	Max_Min_Sum_Char(ISEdge,lnnew,&max_gray,&k,&min_gray);
	average_total=k/lnnew;
	DMA_fieldToRectDMA4ABSync( ISD, ws, 1,0,0, hsn, ws, ISEdge);//hs must be multiple of 8
	DMA_Wait_total();
		//-----------finding threshold for black object detection in low contrast background--------------
	if ( (average_total - min_gray ) > 12 )	// noise level
		thd_b_f = (float) 0.65 * average_total + (float) 0.35 * min_gray-3 ;		// 0.55 => alpha , 0.45 => beta , 5 => gama
	else
		thd_b_f = 0;
	thd_b_old=thd_b_f;
	thd_b = (unsigned char) thd_b_f;
	thd_b=_min2(thd_b,200);
	if ( (max_gray - average_total) > 12 )	// noise level
		thd_w_f = (float) 0.65 * average_total + (float) 0.35 * max_gray ;	// 0.55 => alpha , 0.45 => beta , 5 => gama
	else
		thd_w_f = 255;
	thd_w_old=thd_w_f;
	thd_w = (unsigned char) thd_w_f;
	thd_w=_max2(thd_w,50);
	th_b_3210= (thd_b) | (thd_b<<8) | (thd_b<<16) | (thd_b<<24);
	th_w_3210= (thd_w) | (thd_w<<8) | (thd_w<<16) | (thd_w<<24);
	if (idn_white==1)
	{
		Binary_Func_Bit(ISEdge,ln,th_w_3210,ISEdge);
	}
	else// if (idn_black==1 )*/
	{
		Binary_FuncB_Bit(ISEdge,ln,th_b_3210,ISEdge);
	}
	hsd=hsn/8+1;
	ld=hsd*ws;
	DMA_fieldToRectDMA4ABSync( ISEdge, ws, 1,0,0, hsd, ws, ISD);//hs must be multiple of 8
	DMA_Wait_total();
	free(ISEdge);
	ISEdge=(unsigned char*)malloc(LengthISn/2*sizeof(unsigned char));
	ISEdge2=(unsigned char*)malloc(ld*sizeof(unsigned char));
	ISEdge3=(unsigned char*)malloc(ld*sizeof(unsigned char));
	DMA_fieldToRectDMA4ABSync( ISD, ws, 1,0,0, hsd, ws, ISEdge);//hs must be multiple of 8
	DMA_Wait_total();
	EQ_Func_Mat(ISEdge,ld,ISEdge2);
	Dilate_Bit_Func_total(ISEdge2,m_morphology,ISEdge+ws/8+1,hsn-4,ws);
	Erosion_Bit_Func_total(ISEdge,d_morphology,ISEdge2+ws/8+1,hsn-4,ws);
	m_morphology[0] = mask2[0] > 0 ? ~0U : 0;
	m_morphology[1] = mask2[1] > 0 ? ~0U : 0;
	m_morphology[2] = mask2[2] > 0 ? ~0U : 0;
	m_morphology[3] = mask2[3] > 0 ? ~0U : 0;
	m_morphology[4] = mask2[4] > 0 ? ~0U : 0;
	m_morphology[5] = mask2[5] > 0 ? ~0U : 0;
	m_morphology[6] = mask2[6] > 0 ? ~0U : 0;
	m_morphology[7] = mask2[7] > 0 ? ~0U : 0;
	m_morphology[8] = mask2[8] > 0 ? ~0U : 0;
	Dilate_Bit_Func_total(ISEdge2,m_morphology,ISEdge+ws/8+1,hsn-4,ws);
	Dilate_Bit_Func_total(ISEdge,m_morphology,ISEdge2+ws/8+1,hsn-4,ws);
	for (i=0;i<hsn1;i++)
	{
		ISEdge2[i*ws/8]=0x00;
		ISEdge2[(i+1)*ws/8-1]=0x00;
	}
	for (i=0;i<(ws/8*4);i++)
		ISEdge2[i]=0;
	for (i=(hsn1-4)*ws/8;i<hsn1*ws/8;i++)
		ISEdge2[i]=0;
	DMA_fieldToRectDMA4ABSync( ISEdge2, ws/8, 1,0,0, hsn, ws/8, ISD3);//hs must be multiple of 8
	DMA_Wait_total();
	Resizehalfcol_Bit_Func1_MD(ISEdge2,ISEdge2,hsn,ws);
	ResizeQuart_Bit_Func_MD(ISEdge2,ISEdge2,hsn/2,ws);//ws must be multiple of 32
	//Resizehalfcol_Bit_Func2_MD(ISEdge2,ISEdge2,hsn/2,ws/4);
	Unpack_Bit_Func(ISEdge2,((ws/4)*(hsn/2)/8),ISEdge);
	DMA_fieldToRectDMA4ABSync( ISEdge, ws/4, 1,0,0, hsn1/2, ws2/4, ISD3Resize);//hs must be multiple of 8
	DMA_Wait_total();
	ErosionBitFuncTotal_MD(ISEdge2,d_morphology,ISEdge,hsn1/2,ws/4);
	//Erosion_Bit_Func_total(ISEdge2,d_morphology,ISEdge+ws/8+1,hsn-4,ws);
	NOT_Func(ISEdge,(ws/4)*(hsn/2)/8,ISEdge);
	AND_Func(ISEdge2,ISEdge,ISEdge2,(hsn/2)*(ws/4)/8);
	EQ_Func_Mat(ISEdge2,ld,ISD);
	free(ISEdge);
	free(ISEdge2);
	free(ISEdge3);
	ISEdge=(unsigned char*)malloc(LengthISn*sizeof(char));
	DMA_fieldToRectDMA4ABSync( ISD3, ws, 1,0,0, hsd, ws, IS1);//hs must be multiple of 8
	DMA_Wait_total();
	Unpack_Bit_Func(IS1,(ws*hsn/8),ISEdge);
	DMA_fieldToRectDMA4ABSync( ISEdge, ws, 1,0,0, hsn1, ws2, ISD3);//hs must be multiple of 8
	DMA_Wait_total();
	DMA_fieldToRectDMA4ABSync( ISD, ws, 1,0,0, hsd, ws, IS1);//hs must be multiple of 8
	DMA_Wait_total();
	Unpack_Bit_Func(IS1,(ws*hsn/8),ISEdge);
	DMA_fieldToRectDMA4ABSync( ISEdge, ws/4, 1,0,0, hsn1/2, ws2/4, ISDResize);//hs must be multiple of 8
	DMA_Wait_total();
	free(ISEdge);

}
