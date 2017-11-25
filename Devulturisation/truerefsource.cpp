//
//#include <vector>
//#include <iostream>
//#include <iomanip>
//#include <stdio.h>
//
//using namespace std;
//extern "C" {
//#include "jpeglib.h"
//}
//
//#include <setjmp.h>
//typedef unsigned char u8;
//typedef short i16;
//typedef unsigned short u16;
//
//#define getjsample(value)  ((int) (value))
//#define centerjsample	128
//
//#define dctsize		    8
//#define dctsize2	    64
//
//static const double aanscalefactor[dctsize] =
//{
//	1.0, 1.3870398453221475, 1.3065629648763766, 1.1758756024193588,
//	1.0, 0.7856949583871022, 0.5411961001461971, 0.2758993792829431,
//};
//
//void jpeg_init_fdct_table(const u16 *quantptr, float quant_table[dctsize2])
//{
//	for (int row = 0, i = 0; row < dctsize; row++)
//		for (int col = 0; col < dctsize; i++, col++)
//			quant_table[i] = 0.125 / (aanscalefactor[row] * aanscalefactor[col] * quantptr[i]); // 1/64
//}
//
//void jpeg_init_idct_table(const u16 *quantptr, float quant_table[dctsize2])
//{
//	for (int row = 0, i = 0; row < dctsize; row++)
//		for (int col = 0; col < dctsize; i++, col++)
//			quant_table[i] = quantptr[i] * aanscalefactor[row] * aanscalefactor[col] * 0.125;
//}
//
//u8 range_limit(int x)
//{
//	return x < 0 ? 0 : x > 0xFF ? 0xFF : x;
//}
//
//void jpeg_idct_float(const float *inptr, u8 *output_buf, int output_stride)
//{
//	float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
//	float tmp10, tmp11, tmp12, tmp13;
//	float z5, z10, z11, z12, z13;
//	float * wsptr;
//	int ctr;
//	float workspace[dctsize2];
//
//	wsptr = workspace;
//	for (ctr = dctsize; ctr > 0; ctr--) {
//		if (inptr[dctsize * 1] == 0 && inptr[dctsize * 2] == 0 &&
//			inptr[dctsize * 3] == 0 && inptr[dctsize * 4] == 0 &&
//			inptr[dctsize * 5] == 0 && inptr[dctsize * 6] == 0 &&
//			inptr[dctsize * 7] == 0)
//		{
//			float dcval = inptr[dctsize * 0];
//
//			wsptr[dctsize * 0] = dcval;
//			wsptr[dctsize * 1] = dcval;
//			wsptr[dctsize * 2] = dcval;
//			wsptr[dctsize * 3] = dcval;
//			wsptr[dctsize * 4] = dcval;
//			wsptr[dctsize * 5] = dcval;
//			wsptr[dctsize * 6] = dcval;
//			wsptr[dctsize * 7] = dcval;
//
//			inptr++;
//			wsptr++;
//			continue;
//		}
//
//		tmp0 = inptr[dctsize * 0];
//		tmp1 = inptr[dctsize * 2];
//		tmp2 = inptr[dctsize * 4];
//		tmp3 = inptr[dctsize * 6];
//
//		tmp10 = tmp0 + tmp2;
//		tmp11 = tmp0 - tmp2;
//
//		tmp13 = tmp1 + tmp3;
//		tmp12 = (tmp1 - tmp3) * 1.414213562f - tmp13;
//
//		tmp0 = tmp10 + tmp13;
//		tmp3 = tmp10 - tmp13;
//		tmp1 = tmp11 + tmp12;
//		tmp2 = tmp11 - tmp12;
//
//		tmp4 = inptr[dctsize * 1];
//		tmp5 = inptr[dctsize * 3];
//		tmp6 = inptr[dctsize * 5];
//		tmp7 = inptr[dctsize * 7];
//
//		z13 = tmp6 + tmp5;
//		z10 = tmp6 - tmp5;
//		z11 = tmp4 + tmp7;
//		z12 = tmp4 - tmp7;
//
//		tmp7 = z11 + z13;
//		tmp11 = (z11 - z13) * 1.414213562f;
//
//		z5 = (z10 + z12) * 1.847759065f;
//		tmp10 = z5 - z12 * 1.082392200f;
//		tmp12 = z5 - z10 * 2.613125930f;
//
//		tmp6 = tmp12 - tmp7;
//		tmp5 = tmp11 - tmp6;
//		tmp4 = tmp10 - tmp5;
//
//		wsptr[dctsize * 0] = tmp0 + tmp7;
//		wsptr[dctsize * 7] = tmp0 - tmp7;
//		wsptr[dctsize * 1] = tmp1 + tmp6;
//		wsptr[dctsize * 6] = tmp1 - tmp6;
//		wsptr[dctsize * 2] = tmp2 + tmp5;
//		wsptr[dctsize * 5] = tmp2 - tmp5;
//		wsptr[dctsize * 3] = tmp3 + tmp4;
//		wsptr[dctsize * 4] = tmp3 - tmp4;
//
//		inptr++;
//		wsptr++;
//	}
//
//	/* pass 2: process rows. */
//
//	wsptr = workspace;
//	for (ctr = 0; ctr < dctsize; ctr++) {
//		/* prepare range-limit and float->int conversion */
//		z5 = wsptr[0] + (centerjsample + 0.5f);
//		tmp10 = z5 + wsptr[4];
//		tmp11 = z5 - wsptr[4];
//
//		tmp13 = wsptr[2] + wsptr[6];
//		tmp12 = (wsptr[2] - wsptr[6]) * 1.414213562f - tmp13;
//
//		tmp0 = tmp10 + tmp13;
//		tmp3 = tmp10 - tmp13;
//		tmp1 = tmp11 + tmp12;
//		tmp2 = tmp11 - tmp12;
//
//		z13 = wsptr[5] + wsptr[3];
//		z10 = wsptr[5] - wsptr[3];
//		z11 = wsptr[1] + wsptr[7];
//		z12 = wsptr[1] - wsptr[7];
//
//		tmp7 = z11 + z13;
//		tmp11 = (z11 - z13) * 1.414213562f;
//
//		z5 = (z10 + z12) * 1.847759065f;
//		tmp10 = z5 - z12 * 1.082392200f;
//		tmp12 = z5 - z10 * 2.613125930f;
//
//		tmp6 = tmp12 - tmp7;
//		tmp5 = tmp11 - tmp6;
//		tmp4 = tmp10 - tmp5;
//
//		/* final output stage: float->int conversion and range-limit */
//		output_buf[0] = range_limit((int)(tmp0 + tmp7));
//		output_buf[7] = range_limit((int)(tmp0 - tmp7));
//		output_buf[1] = range_limit((int)(tmp1 + tmp6));
//		output_buf[6] = range_limit((int)(tmp1 - tmp6));
//		output_buf[2] = range_limit((int)(tmp2 + tmp5));
//		output_buf[5] = range_limit((int)(tmp2 - tmp5));
//		output_buf[3] = range_limit((int)(tmp3 + tmp4));
//		output_buf[4] = range_limit((int)(tmp3 - tmp4));
//
//		wsptr += dctsize;
//		output_buf += output_stride;
//	}
//}
//
//
//void jpeg_fdct_float(float *outptr, const u8 *input_buf, int input_stride)
//{
//	float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
//	float tmp10, tmp11, tmp12, tmp13;
//	float z1, z2, z3, z4, z5, z11, z13;
//	float *dataptr;
//	int ctr;
//
//	/* pass 1: process rows. */
//
//	dataptr = outptr;
//	for (ctr = 0; ctr < dctsize; ctr++) {
//		tmp0 = (float)(getjsample(input_buf[0]) + getjsample(input_buf[7]));
//		tmp7 = (float)(getjsample(input_buf[0]) - getjsample(input_buf[7]));
//		tmp1 = (float)(getjsample(input_buf[1]) + getjsample(input_buf[6]));
//		tmp6 = (float)(getjsample(input_buf[1]) - getjsample(input_buf[6]));
//		tmp2 = (float)(getjsample(input_buf[2]) + getjsample(input_buf[5]));
//		tmp5 = (float)(getjsample(input_buf[2]) - getjsample(input_buf[5]));
//		tmp3 = (float)(getjsample(input_buf[3]) + getjsample(input_buf[4]));
//		tmp4 = (float)(getjsample(input_buf[3]) - getjsample(input_buf[4]));
//
//		tmp10 = tmp0 + tmp3;
//		tmp13 = tmp0 - tmp3;
//		tmp11 = tmp1 + tmp2;
//		tmp12 = tmp1 - tmp2;
//
//		/* apply unsigned->signed conversion. */
//		dataptr[0] = tmp10 + tmp11 - 8 * centerjsample;
//		dataptr[4] = tmp10 - tmp11;
//
//		z1 = (tmp12 + tmp13) * 0.707106781f;
//		dataptr[2] = tmp13 + z1;
//		dataptr[6] = tmp13 - z1;
//
//		tmp10 = tmp4 + tmp5;
//		tmp11 = tmp5 + tmp6;
//		tmp12 = tmp6 + tmp7;
//
//		z5 = (tmp10 - tmp12) * 0.382683433f;
//		z2 = 0.541196100f * tmp10 + z5;
//		z4 = 1.306562965f * tmp12 + z5;
//		z3 = tmp11 * 0.707106781f;
//
//		z11 = tmp7 + z3;
//		z13 = tmp7 - z3;
//
//		dataptr[5] = z13 + z2;
//		dataptr[3] = z13 - z2;
//		dataptr[1] = z11 + z4;
//		dataptr[7] = z11 - z4;
//
//		dataptr += dctsize;
//		input_buf += input_stride;
//	}
//
//	/* pass 2: process columns. */
//
//	dataptr = outptr;
//	for (ctr = dctsize - 1; ctr >= 0; ctr--) {
//		tmp0 = dataptr[dctsize * 0] + dataptr[dctsize * 7];
//		tmp7 = dataptr[dctsize * 0] - dataptr[dctsize * 7];
//		tmp1 = dataptr[dctsize * 1] + dataptr[dctsize * 6];
//		tmp6 = dataptr[dctsize * 1] - dataptr[dctsize * 6];
//		tmp2 = dataptr[dctsize * 2] + dataptr[dctsize * 5];
//		tmp5 = dataptr[dctsize * 2] - dataptr[dctsize * 5];
//		tmp3 = dataptr[dctsize * 3] + dataptr[dctsize * 4];
//		tmp4 = dataptr[dctsize * 3] - dataptr[dctsize * 4];
//
//		tmp10 = tmp0 + tmp3;
//		tmp13 = tmp0 - tmp3;
//		tmp11 = tmp1 + tmp2;
//		tmp12 = tmp1 - tmp2;
//
//		dataptr[dctsize * 0] = tmp10 + tmp11;
//		dataptr[dctsize * 4] = tmp10 - tmp11;
//
//		z1 = (tmp12 + tmp13) * 0.707106781f;
//		dataptr[dctsize * 2] = tmp13 + z1;
//		dataptr[dctsize * 6] = tmp13 - z1;
//
//		tmp10 = tmp4 + tmp5;
//		tmp11 = tmp5 + tmp6;
//		tmp12 = tmp6 + tmp7;
//
//		z5 = (tmp10 - tmp12) * 0.382683433f;
//		z2 = 0.541196100f * tmp10 + z5;
//		z4 = 1.306562965f * tmp12 + z5;
//		z3 = tmp11 * 0.707106781f;
//
//		z11 = tmp7 + z3;
//		z13 = tmp7 - z3;
//
//		dataptr[dctsize * 5] = z13 + z2;
//		dataptr[dctsize * 3] = z13 - z2;
//		dataptr[dctsize * 1] = z11 + z4;
//		dataptr[dctsize * 7] = z11 - z4;
//
//		dataptr++;
//	}
//}
//
//char * put_scanline_someplace(JSAMPLE* buffer, int row_stride)
//{
//	char * a = new char[row_stride];
//	for (int i = 0; i < row_stride; i++)
//	{
//		a[i] = buffer[i];
//	}
//	return a;
//
//}
//
//JSAMPLE * image_buffer;
//int image_height;
//int image_width;
//
//
//GLOBAL(void)
//write_JPEG_file(char * filename, int quality)
//{
//
//	struct jpeg_compress_struct cinfo;
//
//	struct jpeg_error_mgr jerr;
//
//	FILE * outfile;
//	JSAMPROW row_pointer[1];
//	int row_stride;
//	cinfo.err = jpeg_std_error(&jerr);
//	jpeg_create_compress(&cinfo);
//
//	if ((outfile = fopen(filename, "wb")) == NULL) {
//		fprintf(stderr, "can't open %s\n", filename);
//		return;
//	}
//	jpeg_stdio_dest(&cinfo, outfile);
//	cinfo.image_width = image_width;
//	cinfo.image_height = image_height;
//	cinfo.input_components = 3;
//	cinfo.in_color_space = JCS_RGB;
//	jpeg_set_defaults(&cinfo);
//
//	jpeg_set_quality(&cinfo, quality, TRUE);
//
//	jpeg_start_compress(&cinfo, TRUE);
//
//	row_stride = image_width * 3;
//
//	while (cinfo.next_scanline < cinfo.image_height) {
//		row_pointer[0] = &image_buffer[cinfo.next_scanline * row_stride];
//		(void)jpeg_write_scanlines(&cinfo, row_pointer, 1);
//	}
//
//
//	jpeg_finish_compress(&cinfo);
//
//	fclose(outfile);
//
//	jpeg_destroy_compress(&cinfo);
//
//}
//
//
//
//struct my_error_mgr {
//	struct jpeg_error_mgr pub;
//
//	jmp_buf setjmp_buffer;
//};
//
//typedef struct my_error_mgr * my_error_ptr;
//
//METHODDEF(void)
//my_error_exit(j_common_ptr cinfo)
//{
//	my_error_ptr myerr = (my_error_ptr)cinfo->err;
//	(*cinfo->err->output_message) (cinfo);
//	longjmp(myerr->setjmp_buffer, 1);
//}
//
//
//void PGMwriting(JSAMPARRAY pich, int w, int h, char *filename, int num) 
//{
//  FILE * infile;
//  int n = strlen(filename);
//  char *newname = new char[n+1];
//  sprintf(newname, "%.*s%i.pgm",n-4,filename, num);
//  infile = fopen(newname, "wb");
//  fprintf(infile, "P5\n%i %i\n255\n", w, h);
//  for (int i = 0; i < (int)(h); i++)
//  {
//    fwrite(pich[i], 1, w, infile);
//  }
//  fclose(infile);
//}
//
//
//
//
//GLOBAL(int)
//read_JPEG_file(char * filename, u8 *pich, int w, int h, float fquant[dctsize2], float iquant[dctsize2])
//{
//	struct jpeg_decompress_struct cinfo;
//	struct my_error_mgr jerr;
//	FILE * infile;
//  JSAMPARRAY buffer;
//	int row_stride;
//
//	if ((infile = fopen(filename, "rb")) == NULL) {
//		fprintf(stderr, "can't open %s\n", filename);
//		return 0;
//	}
//	cinfo.err = jpeg_std_error(&jerr.pub);
//	jerr.pub.error_exit = my_error_exit;
//	if (setjmp(jerr.setjmp_buffer)) {
//		jpeg_destroy_decompress(&cinfo);
//		fclose(infile);
//		return 0;
//	}
//
//	jpeg_create_decompress(&cinfo);
//	jpeg_stdio_src(&cinfo, infile);
//	(void)jpeg_read_header(&cinfo, TRUE);
//	cinfo.out_color_space = JCS_GRAYSCALE;
//  //cinfo.raw_data_out = TRUE;
//	(void)jpeg_start_decompress(&cinfo);
//	row_stride = cinfo.output_width * cinfo.output_components;
//
//	buffer = (*cinfo.mem->alloc_sarray)
//		((j_common_ptr)&cinfo, JPOOL_IMAGE, row_stride, 1);
//
//
//	u16 quant[dctsize2];
//	for (int i = 0; i < dctsize2; i++) {
//		quant[i] = 1;
//	}
//	///uint16 -> u16
//	JQUANT_TBL** quants = cinfo.quant_tbl_ptrs;
//	for (int i = 0; i < dctsize2; i++) {
//			quant[i] = quants[0]->quantval[i];
//			cout << setw(4) << quants[0]->quantval[i]<<' ';
//      if ((i + 1) % 8 == 0) {
//        cout << endl;
//      }
//			//quant[i] = (u16)quants[i];
//	}
//
//	//float fquant[dctsize2], iquant[dctsize2];
//	jpeg_init_fdct_table(quant, fquant);
//	jpeg_init_idct_table(quant, iquant);
//
//	u8 buf1[dctsize2], buf2[dctsize2];
//  printf("%i\n", cinfo.comp_info->height_in_blocks*dctsize);
//  printf("%i\n", cinfo.comp_info->width_in_blocks*dctsize);
//  printf("%i x %i \n", cinfo.max_h_samp_factor, cinfo.max_v_samp_factor);
//  printf("%i x %i \n", cinfo.min_DCT_h_scaled_size, cinfo.min_DCT_v_scaled_size);
//  printf("%i x %i \n", cinfo.output_height, cinfo.output_width);
//  w = cinfo.output_width;
//  h = cinfo.output_height;
//	JSAMPARRAY a = new JSAMPROW[cinfo.output_height];
//	//JSAMPARRAY resDCT = new JSAMPROW[cinfo.output_height];
//	//int** preResDCT = new int*[cinfo.output_height];
//	//int **b= new int* [cinfo.output_height];
//  pich = new u8[w*h*cinfo.output_components];
//	int i = 0;
//	while (cinfo.output_scanline < cinfo.output_height) {
//		(void)jpeg_read_scanlines(&cinfo, buffer, 1);
//		a[i] = new JSAMPLE[row_stride];
//		//resDCT[i] = new JSAMPLE[row_stride];
//		//preResDCT[i]= new int[row_stride];
//		//b[i] = new int[row_stride];
//		for (int j = 0; j < row_stride; j++)
//		{
//			a[i][j] = buffer[0][j];
//			//preResDCT[i][j] = 0;
//			//resDCT[i][j] = 0;
//			//b[i][j] = 0;
//		}
//		i++;
//	}
//  for (int y = 0; y < cinfo.output_height; y++)
//  {
//    for (int z = 0; z < row_stride; z++)
//    {
//      pich[y*row_stride + z] = a[y][z];
//    }
//  }
//
//  PGMwriting(a, cinfo.output_width, cinfo.output_height, filename, 0);
//	//float  tmp[64] ;
//	//unsigned char  tmp2[64];
//	//for (i = 0; i < (int)(cinfo.output_height)-7; i++)
//	//{
//	//	for (int j = 0; j < (int)(row_stride)-7; j++)
//	//	{
//	//		
//	//		for (int y = 0; y < 8; y++)
//	//		{
//	//			for (int z = 0; z < 8; z++)
//	//			{
//	//				tmp2[y * 8 + z] = a[i+y][j+z];
//	//				
//	//			}
//	//		}
//	//		jpeg_fdct_float(tmp, tmp2, dctsize);
//	//		for (int i = 0; i < dctsize2; i++) {
//	//			tmp[i] = round(tmp[i] * fquant[i]) * iquant[i];
//	//		}
//
//	//		jpeg_idct_float(tmp, tmp2, dctsize);
//	//		for (int y = 0; y < 8; y++)
//	//		{
//	//			for (int z = 0; z < 8; z++)
//	//			{
//	//				preResDCT[y + i][z + j]+= tmp2[y * 8 + z];
//	//				b[y + i][z + j]++;
//	//			}
//	//		}
//	//	}
//	//}
//
//
//	//for (int y = 0; y < cinfo.output_height; y++)
//	//{
//	//	for (int z = 0; z < row_stride; z++)
//	//	{
//	//		
//	//		preResDCT[y][z] /= b[y][z];
//	//		resDCT[y][z] = range_limit(preResDCT[y][z]);
//	//	}
//	//}
// // PGMwriting(resDCT, cinfo.output_width, cinfo.output_height, filename, 1);
//
//	(void)jpeg_finish_decompress(&cinfo);
//	jpeg_destroy_decompress(&cinfo);
//	
//	return 1;
//}
