
#include <vector>
#include <iostream>
using namespace std;
#include <stdio.h>


extern "C" {

#include "jpeglib.h"

}

static const double aanscalefactor[DCTSIZE]= {
	1.0, 1.387039845, 1.306562965, 1.175875602,
	1.0, 0.785694958, 0.541196100, 0.275899379
};
double superS[] = { 64, 88.7706, 83.6201, 75.2559, 64, 50.2844, 34.6366, 17.6576,
88.7705, 123.128 ,115.984 ,104.383, 88.7707, 69.7467 ,48.0424, 24.4917,
83.6199, 115.984, 109.255, 98.3268, 83.62, 65.6998, 45.2548, 23.0707,
75.256, 104.383, 98.3268, 88.4918, 75.2558, 59.1281, 40.7283, 20.763,
63.9998, 88.7706, 83.62, 75.2561, 64, 50.2845, 34.6366, 17.6575,
50.2845, 69.7467, 65.6998, 59.1282, 50.2845, 39.5083, 27.2138, 13.8735,
34.6365, 48.0423, 45.2548, 40.7282, 34.6366, 27.2138 ,18.7452, 9.55619,
17.6576, 24.4917, 23.0707, 20.7631, 17.6575 ,13.8735, 9.5562, 4.87171 };
#include <setjmp.h>
//void createSuperS()
//{
//	for (int i = 0; i < 8; i++)
//	{
//		for (int j = 0; j < 8; j++)
//		{
//			superS[i][j] = aanscalefactor[j]* aanscalefactor[i];
//		}
//	}
//}
void jpeg_fdct_float(float * data)
{
	float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
	float tmp10, tmp11, tmp12, tmp13;
	float z1, z2, z3, z4, z5, z11, z13;
	float *dataptr;
	int ctr;

	/* Pass 1: process rows. */

	dataptr = data;
	for (ctr = DCTSIZE - 1; ctr >= 0; ctr--) {
		tmp0 = dataptr[0] + dataptr[7];
		tmp7 = dataptr[0] - dataptr[7];
		tmp1 = dataptr[1] + dataptr[6];
		tmp6 = dataptr[1] - dataptr[6];
		tmp2 = dataptr[2] + dataptr[5];
		tmp5 = dataptr[2] - dataptr[5];
		tmp3 = dataptr[3] + dataptr[4];
		tmp4 = dataptr[3] - dataptr[4];

		/* Even part */

		tmp10 = tmp0 + tmp3;	/* phase 2 */
		tmp13 = tmp0 - tmp3;
		tmp11 = tmp1 + tmp2;
		tmp12 = tmp1 - tmp2;


		dataptr[0] = tmp10 + tmp11; /* phase 3 */
		dataptr[4] = tmp10 - tmp11;

		z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
		dataptr[2] =( tmp13 + z1);	/* phase 5 */
		dataptr[6] =( tmp13 - z1);

		/* Odd part */

		tmp10 = tmp4 + tmp5;	/* phase 2 */
		tmp11 = tmp5 + tmp6;
		tmp12 = tmp6 + tmp7;

		/* The rotator is modified from fig 4-8 to avoid extra negations. */
		z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
		z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
		z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
		z3 = tmp11 * ((float) 0.707106781); /* c4 */

		z11 = tmp7 + z3;		/* phase 5 */
		z13 = tmp7 - z3;

		dataptr[5] = (z13 + z2);	/* phase 6 */
		dataptr[3] =(z13 - z2);
		dataptr[1] = (z11 + z4);
		dataptr[7] = (z11 - z4);

		dataptr += DCTSIZE;		/* advance pointer to next row */
	}

	/* Pass 2: process columns. */

	dataptr = data;
	for (ctr = DCTSIZE - 1; ctr >= 0; ctr--) {
		tmp0 = dataptr[DCTSIZE * 0] + dataptr[DCTSIZE * 7];
		tmp7 = dataptr[DCTSIZE * 0] - dataptr[DCTSIZE * 7];
		tmp1 = dataptr[DCTSIZE * 1] + dataptr[DCTSIZE * 6];
		tmp6 = dataptr[DCTSIZE * 1] - dataptr[DCTSIZE * 6];
		tmp2 = dataptr[DCTSIZE * 2] + dataptr[DCTSIZE * 5];
		tmp5 = dataptr[DCTSIZE * 2] - dataptr[DCTSIZE * 5];
		tmp3 = dataptr[DCTSIZE * 3] + dataptr[DCTSIZE * 4];
		tmp4 = dataptr[DCTSIZE * 3] - dataptr[DCTSIZE * 4];

		/* Even part */

		tmp10 = tmp0 + tmp3;	/* phase 2 */
		tmp13 = tmp0 - tmp3;
		tmp11 = tmp1 + tmp2;
		tmp12 = tmp1 - tmp2;

		dataptr[DCTSIZE * 0] = tmp10 + tmp11; /* phase 3 */
		dataptr[DCTSIZE * 4] = tmp10 - tmp11;

		z1 = (tmp12 + tmp13) * ((float) 0.707106781); /* c4 */
		dataptr[DCTSIZE * 2] = tmp13 + z1; /* phase 5 */
		dataptr[DCTSIZE * 6] = tmp13 - z1;

		/* Odd part */

		tmp10 = tmp4 + tmp5;	/* phase 2 */
		tmp11 = tmp5 + tmp6;
		tmp12 = tmp6 + tmp7;

		/* The rotator is modified from fig 4-8 to avoid extra negations. */
		z5 = (tmp10 - tmp12) * ((float) 0.382683433); /* c6 */
		z2 = ((float) 0.541196100) * tmp10 + z5; /* c2-c6 */
		z4 = ((float) 1.306562965) * tmp12 + z5; /* c2+c6 */
		z3 = tmp11 * ((float) 0.707106781); /* c4 */

		z11 = tmp7 + z3;		/* phase 5 */
		z13 = tmp7 - z3;

		dataptr[DCTSIZE * 5] = z13 + z2; /* phase 6 */
		dataptr[DCTSIZE * 3] = z13 - z2;
		dataptr[DCTSIZE * 1] = z11 + z4;
		dataptr[DCTSIZE * 7] = z11 - z4;





	
		dataptr++;			/* advance pointer to next column */
	}
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			data[i*8+j] /=superS[i*8+j];
		}
	}
	
	
}

GLOBAL(void)
jpeg_idct_float(float * in)
{
	float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
	float tmp10, tmp11, tmp12, tmp13;
	float z5, z10, z11, z12, z13;
	float* inptr;

	
	int ctr;
	float workspace[DCTSIZE2]; /* buffers data between passes */


		/* Pass 1: process columns from input, store into work array. */

		inptr = in;

	for (ctr = DCTSIZE; ctr > 0; ctr--) {
		/* Due to quantization, we will usually find that many of the input
		* coefficients are zero, especially the AC terms.  We can exploit this
		* by short-circuiting the IDCT calculation for any column in which all
		* the AC terms are zero.  In that case each output is equal to the
		* DC coefficient (with scale factor as needed).
		* With typical images and quantization tables, half or more of the
		* column DCT calculations can be simplified this way.
		*/

		/* Even part */

		tmp0 = inptr[DCTSIZE * 0];
		tmp1 = inptr[DCTSIZE * 2];
		tmp2 = inptr[DCTSIZE * 4];
		tmp3 = inptr[DCTSIZE * 6];

		tmp10 = tmp0 + tmp2;	/* phase 3 */
		tmp11 = tmp0 - tmp2;

		tmp13 = tmp1 + tmp3;	/* phases 5-3 */
		tmp12 = (tmp1 - tmp3) * ((float) 1.414213562) - tmp13; /* 2*c4 */

		tmp0 = tmp10 + tmp13;	/* phase 2 */
		tmp3 = tmp10 - tmp13;
		tmp1 = tmp11 + tmp12;
		tmp2 = tmp11 - tmp12;

		/* Odd part */

		tmp4 = inptr[DCTSIZE * 1];
		tmp5 = inptr[DCTSIZE * 3];
		tmp6 = inptr[DCTSIZE * 5];
		tmp7 = inptr[DCTSIZE * 7];

		z13 = tmp6 + tmp5;		/* phase 6 */
		z10 = tmp6 - tmp5;
		z11 = tmp4 + tmp7;
		z12 = tmp4 - tmp7;

		tmp7 = z11 + z13;		/* phase 5 */
		tmp11 = (z11 - z13) * ((float) 1.414213562); /* 2*c4 */

		z5 = (z10 + z12) * ((float) 1.847759065); /* 2*c2 */
		tmp10 = ((float) 1.082392200) * z12 - z5; /* 2*(c2-c6) */
		tmp12 = ((float)-2.613125930) * z10 + z5; /* -2*(c2+c6) */

		tmp6 = tmp12 - tmp7;	/* phase 2 */
		tmp5 = tmp11 - tmp6;
		tmp4 = tmp10 + tmp5;

		inptr[DCTSIZE * 0] = tmp0 + tmp7;
		inptr[DCTSIZE * 7] = tmp0 - tmp7;
		inptr[DCTSIZE * 1] = tmp1 + tmp6;
		inptr[DCTSIZE * 6] = tmp1 - tmp6;
		inptr[DCTSIZE * 2] = tmp2 + tmp5;
		inptr[DCTSIZE * 5] = tmp2 - tmp5;
		inptr[DCTSIZE * 4] = tmp3 + tmp4;
		inptr[DCTSIZE * 3] = tmp3 - tmp4;

		inptr++;			/* advance pointers to next column */
		
		
	}

	/* Pass 2: process rows from work array, store into output array. */
	/* Note that we must descale the results by a factor of 8 == 2**3. */

	inptr = workspace;
	for (ctr = 0; ctr < DCTSIZE; ctr++) {
		/* Rows of zeroes can be exploited in the same way as we did with columns.
		* However, the column calculation has created many nonzero AC terms, so
		* the simplification applies less often (typically 5% to 10% of the time).
		* And testing floats for zero is relatively expensive, so we don't bother.
		*/

		/* Even part */

		tmp10 = inptr[0] + inptr[4];
		tmp11 = inptr[0] - inptr[4];

		tmp13 = inptr[2] + inptr[6];
		tmp12 = (inptr[2] - inptr[6]) * ((float) 1.414213562) - tmp13;

		tmp0 = tmp10 + tmp13;
		tmp3 = tmp10 - tmp13;
		tmp1 = tmp11 + tmp12;
		tmp2 = tmp11 - tmp12;

		/* Odd part */

		z13 = inptr[5] + inptr[3];
		z10 = inptr[5] - inptr[3];
		z11 = inptr[1] + inptr[7];
		z12 = inptr[1] - inptr[7];

		tmp7 = z11 + z13;
		tmp11 = (z11 - z13) * ((float) 1.414213562);

		z5 = (z10 + z12) * ((float) 1.847759065); /* 2*c2 */
		tmp10 = ((float) 1.082392200) * z12 - z5; /* 2*(c2-c6) */
		tmp12 = ((float)-2.613125930) * z10 + z5; /* -2*(c2+c6) */

		tmp6 = tmp12 - tmp7;
		tmp5 = tmp11 - tmp6;
		tmp4 = tmp10 + tmp5;

		/* Final output stage: scale down by a factor of 8 and range-limit */

		inptr[0] = tmp0 + tmp7;
			
		inptr[7] = tmp0 - tmp7;

		inptr[1] = tmp1 + tmp6;

		inptr[6] = tmp1 - tmp6;

		inptr[2] = tmp2 + tmp5;

		inptr[5] = tmp2 - tmp5;

		inptr[4] = tmp3 + tmp4;
		
		inptr[3] = tmp3 - tmp4;
		

		inptr ++;		/* advance pointer to next row */
	}
}

char * put_scanline_someplace(JSAMPLE* buffer, int row_stride)
{
	char * a = new char[row_stride];
	for (int i = 0; i < row_stride; i++)
	{
		a[i] = buffer[i];
	}
	return a;

}

JSAMPLE * image_buffer;
int image_height;
int image_width;


GLOBAL(void)
write_JPEG_file(char * filename, int quality)
{
	
	struct jpeg_compress_struct cinfo;

	struct jpeg_error_mgr jerr;

	FILE * outfile;
	JSAMPROW row_pointer[1];
	int row_stride;	
	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	if ((outfile = fopen(filename, "wb")) == NULL) {
		fprintf(stderr, "can't open %s\n", filename);
		return;
	}
	jpeg_stdio_dest(&cinfo, outfile);
	cinfo.image_width = image_width;
	cinfo.image_height = image_height;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;
	jpeg_set_defaults(&cinfo);

	jpeg_set_quality(&cinfo, quality, TRUE );

	jpeg_start_compress(&cinfo, TRUE);

	row_stride = image_width * 3;

	while (cinfo.next_scanline < cinfo.image_height) {
		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
		(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}


	jpeg_finish_compress(&cinfo);

	fclose(outfile);

	jpeg_destroy_compress(&cinfo);

}




struct my_error_mgr {
	struct jpeg_error_mgr pub;

	jmp_buf setjmp_buffer;
};

typedef struct my_error_mgr * my_error_ptr;



METHODDEF(void)
my_error_exit(j_common_ptr cinfo)
{
	my_error_ptr myerr = (my_error_ptr)cinfo->err;
	(*cinfo->err->output_message) (cinfo);
	longjmp(myerr->setjmp_buffer, 1);
}

GLOBAL(int)
read_JPEG_file(char * filename)
{
	struct jpeg_decompress_struct cinfo;
	struct my_error_mgr jerr;
	FILE * infile;	
	JSAMPARRAY buffer;
	int row_stride;	

	if ((infile = fopen(filename, "rb")) == NULL) {
		fprintf(stderr, "can't open %s\n", filename);
		return 0;
	}
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;
	if (setjmp(jerr.setjmp_buffer)) {
		jpeg_destroy_decompress(&cinfo);
		fclose(infile);
		return 0;
	}
	jpeg_create_decompress(&cinfo);
	jpeg_stdio_src(&cinfo, infile);
	(void)jpeg_read_header(&cinfo, TRUE);
	cinfo.out_color_space = JCS_GRAYSCALE;
	(void)jpeg_start_decompress(&cinfo);
	row_stride = cinfo.output_width * cinfo.output_components;

	buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
	infile = fopen("test.pgm", "wb");
	fprintf(infile, "P5\n%i %i\n255\n", cinfo.output_width, cinfo.output_height);
	
	JSAMPARRAY a = new JSAMPROW[cinfo.output_height];
	
	int i = 0;
	while (cinfo.output_scanline < cinfo.output_height) {
		(void)jpeg_read_scanlines(&cinfo, buffer, 1);
		fwrite(buffer[0], 1, row_stride, infile);
		a[i] = new JSAMPLE[row_stride];
		for (int j = 0; j < row_stride; j++)
		{
			a[i][j] = buffer[0][j];
		}
		i++;
		
	}
	float ab[] = { -16342, 2084, -10049, 10117, 2786, -659, -4905, 12975,
		10579, 8081, -10678, 11762, 6898, 444, -6422, -15892,
		-13388, -4441, -11556, -10947, 16008, -1779, -12481, -16230,
		-16091, -4001, 1038, 2333, 3335, 3512, -10936, 5343,
		-1612, -4845, -14514, 3529, 9284, 9916, 652, -6489,
		12320, 7428, 14939, 13950, 1290, -11719, -1242, -8672,
		11870, -9515, 9164, 11261, 16279, 16374, 3654, -3524,
		-7660, -6642, 11146, -15605, -4067, -13348, 5807, -14541 };

	for ( i = 0; i < 1/*(int)(cinfo.output_height / 8)*/; i++)
	{
		for (int j = 0; j < 1/*(int)(row_stride / 8)*/; j++)
		{
			float * tmp=new float[8*8];
			for (int y = 0; y < 8; y++)
			{
				for (int z = 0; z < 8; z++)
				{
					tmp[y * 8 + z] = (float)(ab[y * 8 + z]);
					cout << tmp[y * 8 + z] << ' ';
				}
				cout << endl;
			}
			cout << endl;
			cout << endl;
			jpeg_fdct_float(tmp);
			jpeg_idct_float(tmp);
			for (int y = 0; y < 8; y++)
			{
				for (int z = 0; z < 8; z++)
				{
					cout << tmp[y * 8 + z]<<' ';
				}
				cout << endl;
			}
		}
	}
	(void)jpeg_finish_decompress(&cinfo);
	jpeg_destroy_decompress(&cinfo);
	fclose(infile);
	return 1;
}



void main() {
	//createSuperS();
/*
	for (int i = 0; i < 8; i++)
	{
		for (int j = 0; j < 8; j++)
		{
			cout<< superS[i][j]<<' ';
		}
		cout << endl;
	}
*/
	read_JPEG_file("H:\\ConsoleApplication5\\jpeg\\testimg.jpg");
}
