#define DCTSIZE		    8
#define DCTSIZE2	    64
#define LOSI          16

#define CENTERJSAMPLE	128

#define range_limit convert_uchar_sat_rtz

void jpeg_idct_float(const float *inptr, float *output_buf, int output_stride)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z5, z10, z11, z12, z13;
  float * wsptr;
  int ctr;
  //float workspace[DCTSIZE2];

  wsptr = output_buf;
  for (ctr = DCTSIZE; ctr > 0; ctr--) {
    tmp0 = inptr[DCTSIZE * 0];
    tmp1 = inptr[DCTSIZE * 2];
    tmp2 = inptr[DCTSIZE * 4];
    tmp3 = inptr[DCTSIZE * 6];

    tmp10 = tmp0 + tmp2;
    tmp11 = tmp0 - tmp2;

    tmp13 = tmp1 + tmp3;
    tmp12 = (tmp1 - tmp3) * 1.414213562f - tmp13;

    tmp0 = tmp10 + tmp13;
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    tmp4 = inptr[DCTSIZE * 1];
    tmp5 = inptr[DCTSIZE * 3];
    tmp6 = inptr[DCTSIZE * 5];
    tmp7 = inptr[DCTSIZE * 7];

    z13 = tmp6 + tmp5;
    z10 = tmp6 - tmp5;
    z11 = tmp4 + tmp7;
    z12 = tmp4 - tmp7;

    tmp7 = z11 + z13;
    tmp11 = (z11 - z13) * 1.414213562f;

    z5 = (z10 + z12) * 1.847759065f;
    tmp10 = z5 - z12 * 1.082392200f;
    tmp12 = z5 - z10 * 2.613125930f;

    tmp6 = tmp12 - tmp7;
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;

    wsptr[DCTSIZE * 0] = tmp0 + tmp7;
    wsptr[DCTSIZE * 7] = tmp0 - tmp7;
    wsptr[DCTSIZE * 1] = tmp1 + tmp6;
    wsptr[DCTSIZE * 6] = tmp1 - tmp6;
    wsptr[DCTSIZE * 2] = tmp2 + tmp5;
    wsptr[DCTSIZE * 5] = tmp2 - tmp5;
    wsptr[DCTSIZE * 3] = tmp3 + tmp4;
    wsptr[DCTSIZE * 4] = tmp3 - tmp4;

    inptr++;
    wsptr++;
  }

  /* pass 2: process rows. */

  wsptr = output_buf;
#pragma unroll
  for (ctr = 0; ctr < DCTSIZE; ctr++) {
    /* prepare range-limit and float->int conversion */
    z5 = wsptr[0] + (CENTERJSAMPLE + 0.5f);
    tmp10 = z5 + wsptr[4];
    tmp11 = z5 - wsptr[4];

    tmp13 = wsptr[2] + wsptr[6];
    tmp12 = (wsptr[2] - wsptr[6]) * 1.414213562f - tmp13;

    tmp0 = tmp10 + tmp13;
    tmp3 = tmp10 - tmp13;
    tmp1 = tmp11 + tmp12;
    tmp2 = tmp11 - tmp12;

    z13 = wsptr[5] + wsptr[3];
    z10 = wsptr[5] - wsptr[3];
    z11 = wsptr[1] + wsptr[7];
    z12 = wsptr[1] - wsptr[7];

    tmp7 = z11 + z13;
    tmp11 = (z11 - z13) * 1.414213562f;

    z5 = (z10 + z12) * 1.847759065f;
    tmp10 = z5 - z12 * 1.082392200f;
    tmp12 = z5 - z10 * 2.613125930f;

    tmp6 = tmp12 - tmp7;
    tmp5 = tmp11 - tmp6;
    tmp4 = tmp10 - tmp5;

    /* final output stage: float->int conversion and range-limit */
    output_buf[0] = range_limit(tmp0 + tmp7);
    output_buf[7] = range_limit(tmp0 - tmp7);
    output_buf[1] = range_limit(tmp1 + tmp6);
    output_buf[6] = range_limit(tmp1 - tmp6);
    output_buf[2] = range_limit(tmp2 + tmp5);
    output_buf[5] = range_limit(tmp2 - tmp5);
    output_buf[3] = range_limit(tmp3 + tmp4);
    output_buf[4] = range_limit(tmp3 - tmp4);

    wsptr += DCTSIZE;
    output_buf += output_stride;
  }
}

void jpeg_fdct_float(float *outptr, __local const float *input_buf, int input_stride)
{
  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z1, z2, z3, z4, z5, z11, z13;
  float *dataptr;
  int ctr;

  /* pass 1: process rows. */

  dataptr = outptr;

#pragma unroll
  for (ctr = 0; ctr < DCTSIZE; ctr++) {
    tmp0 = (input_buf[0] + input_buf[7]);
    tmp7 = (input_buf[0] - input_buf[7]);
    tmp1 = (input_buf[1] + input_buf[6]);
    tmp6 = (input_buf[1] - input_buf[6]);
    tmp2 = (input_buf[2] + input_buf[5]);
    tmp5 = (input_buf[2] - input_buf[5]);
    tmp3 = (input_buf[3] + input_buf[4]);
    tmp4 = (input_buf[3] - input_buf[4]);

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    /* apply unsigned->signed conversion. */
    dataptr[0] = tmp10 + tmp11 - 8 * CENTERJSAMPLE;
    dataptr[4] = tmp10 - tmp11;

    z1 = (tmp12 + tmp13) * 0.707106781f;
    dataptr[2] = tmp13 + z1;
    dataptr[6] = tmp13 - z1;

    tmp10 = tmp4 + tmp5;
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    z5 = (tmp10 - tmp12) * 0.382683433f;
    z2 = 0.541196100f * tmp10 + z5;
    z4 = 1.306562965f * tmp12 + z5;
    z3 = tmp11 * 0.707106781f;

    z11 = tmp7 + z3;
    z13 = tmp7 - z3;

    dataptr[5] = z13 + z2;
    dataptr[3] = z13 - z2;
    dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;

    dataptr += DCTSIZE;
    input_buf += input_stride;
  }

  /* pass 2: process columns. */

  dataptr = outptr;
  for (ctr = DCTSIZE - 1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE * 0] + dataptr[DCTSIZE * 7];
    tmp7 = dataptr[DCTSIZE * 0] - dataptr[DCTSIZE * 7];
    tmp1 = dataptr[DCTSIZE * 1] + dataptr[DCTSIZE * 6];
    tmp6 = dataptr[DCTSIZE * 1] - dataptr[DCTSIZE * 6];
    tmp2 = dataptr[DCTSIZE * 2] + dataptr[DCTSIZE * 5];
    tmp5 = dataptr[DCTSIZE * 2] - dataptr[DCTSIZE * 5];
    tmp3 = dataptr[DCTSIZE * 3] + dataptr[DCTSIZE * 4];
    tmp4 = dataptr[DCTSIZE * 3] - dataptr[DCTSIZE * 4];

    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;

    dataptr[DCTSIZE * 0] = tmp10 + tmp11;
    dataptr[DCTSIZE * 4] = tmp10 - tmp11;

    z1 = (tmp12 + tmp13) * 0.707106781f;
    dataptr[DCTSIZE * 2] = tmp13 + z1;
    dataptr[DCTSIZE * 6] = tmp13 - z1;

    tmp10 = tmp4 + tmp5;
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    z5 = (tmp10 - tmp12) * 0.382683433f;
    z2 = 0.541196100f * tmp10 + z5;
    z4 = 1.306562965f * tmp12 + z5;
    z3 = tmp11 * 0.707106781f;

    z11 = tmp7 + z3;
    z13 = tmp7 - z3;

    dataptr[DCTSIZE * 5] = z13 + z2;
    dataptr[DCTSIZE * 3] = z13 - z2;
    dataptr[DCTSIZE * 1] = z11 + z4;
    dataptr[DCTSIZE * 7] = z11 - z4;

    dataptr++;
  }
}
//nobody knows about my puns there cause of UTF8 -> C1251
__kernel void devultur_most(const __global uchar* pic, __global ushort* res, const __global float *fquant, const __global float *iquant, int rsize, int offset)
{
  uint l64 = get_local_id(0);
  //срез для блока 8*8
  uint lr = l64 / DCTSIZE;
  uint lc = l64 % DCTSIZE;
  //срез для блока 16*4
  uint lrn = l64 / 4;
  uint lcn = l64 % 4;
  //казалось бы, индекс первого элемента внутри данных одной локальной группы
  uint gr = get_group_id(1)*LOSI;
  uint gc = get_group_id(0)*LOSI;
  //вот где-то тут надо задуматься о смещении указателей на определенный шаг

  float ctblock1[DCTSIZE2];

  __local	union {
    float t[LOSI*LOSI];
    float f[LOSI*LOSI];
  }tot;

  //а тут - запись из глобальной ? Без сознания же писалось
  float4 tmptot1 = convert_float4(vload4(0, pic + (gr + lrn)*rsize + gc + lcn * 4 + offset));
  vstore4(tmptot1, 0, &tot.t[lrn*LOSI + lcn * 4]);
  barrier(CLK_LOCAL_MEM_FENCE);

  //2. DCT-QUANT-IDCT
  jpeg_fdct_float(ctblock1, &tot.t[lr*LOSI + lc], LOSI);
#pragma unroll
  for (int i = 0; i < DCTSIZE2; i++) {
    ctblock1[i] = round(ctblock1[i] * fquant[i]) * iquant[i];
  }
  jpeg_idct_float(ctblock1, ctblock1, DCTSIZE);

  barrier(CLK_LOCAL_MEM_FENCE);
  vstore4((float4)(0), 0, &tot.f[lrn*LOSI + lcn * 4]);
  barrier(CLK_LOCAL_MEM_FENCE);

  //3. запись-суммирование с преградами и конями внутри локальной группы
#pragma unroll
  for (int y = 0; y < DCTSIZE; y++)
  {
#pragma unroll
    for (int z = 0; z < DCTSIZE; z++)
    {
      tot.f[(y + lr)*LOSI + z + lc] += ctblock1[y * DCTSIZE + z];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }
  //4. запись данных аккуратно в глоб
  //недопись - допеши
  //ushort4 tmptotfin = vload4(0, res +(gr + lrn)*rsize + gc + lcn * 4 + offset);
  //vstore4(tmptotfin + convert_ushort4_sat_rtz(vload4(0, &fintot[lrn*LOSI+lcn*4])),0, res + (gr + lrn)*rsize + gc + lcn * 4 + offset);

  *(__global ushort4 *)(res + (gr + lrn)*rsize + gc + lcn * 4 + offset) += convert_ushort4_sat_rtz(*(__local float4 *)&tot.f[lrn*LOSI + lcn * 4]);

  //5. ...
  //6. PROFIT!
}

__kernel void devultur_vert(const __global uchar* pic, __global ushort* res, const __global float *fquant, const __global float *iquant, int rsize, int csize, int offset)
{
  uint gr = get_group_id(1)*LOSI;
  uint lr = get_local_id(1);

  float ctblock1[DCTSIZE2];

  __local	union {
    float t[DCTSIZE*LOSI];
    float f[DCTSIZE*LOSI];
  }tot;

  float8 tmptot1 = convert_float8(vload8(0, pic + (gr + lr)*rsize + offset));
  vstore8(tmptot1, 0, &tot.t[lr * 8]);
  tmptot1 = convert_float8(vload8(0, pic + (gr + lr + 8)* rsize + offset));
  vstore8(tmptot1, 0, &tot.t[lr * 8 + DCTSIZE2]);


  jpeg_fdct_float(ctblock1, &tot.t[lr * 8], DCTSIZE);
#pragma unroll
  for (int i = 0; i < DCTSIZE2; i++) {
    ctblock1[i] = round(ctblock1[i] * fquant[i]) * iquant[i];
  }

  jpeg_idct_float(ctblock1, ctblock1, DCTSIZE);

  barrier(CLK_LOCAL_MEM_FENCE);
  vstore16((float16)(0), 0, &tot.f[lr * 16]);
  barrier(CLK_LOCAL_MEM_FENCE);
#pragma unroll
  for (int y = 0; y < DCTSIZE; y++)
  {
#pragma unroll
    for (int z = 0; z < DCTSIZE; z++)
    {
      tot.f[lr * 8 + z + y * 8] += ctblock1[y * 8 + z];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  //ushort8 tmptotfin1 = vload8(0, res + (gr + lr)*rsize + offset) + convert_ushort8_sat_rtz(vload8(0, &tot.f[lr * 8]));
  //vstore8(tmptotfin1, 0, res + (gr + lr)*rsize + offset);
  //tmptotfin1 = vload8(0, res + (gr + lr + 8)* rsize + offset) + convert_ushort8_sat_rtz(vload8(0, &tot.f[lr * 8 + DCTSIZE2]));
  //vstore8(tmptotfin1, 0, res + (gr + lr + 8)* rsize + offset);

  *(__global ushort8 *)(res + (gr + lr)*rsize + offset) += convert_ushort8_sat_rtz(*(__local float8 *)&tot.f[lr * 8]);
  *(__global ushort8 *)(res + (gr + lr + 8)* rsize + offset) += convert_ushort8_sat_rtz(*(__local float8 *)&tot.f[lr * 8 + DCTSIZE2]);
}

__kernel void devultur_hori(const __global uchar* pic, __global ushort* res, const __global float *fquant, const __global float *iquant, int rsize, int csize, int offset)
{
  uint gc = get_group_id(0)*LOSI;
  uint lr = get_local_id(0);

  float ctblock1[DCTSIZE2];

  __local	union {
    float t[DCTSIZE*LOSI];
    float f[DCTSIZE*LOSI];
  }tot;

  float8 tmptot1 = convert_float8(vload8(0, pic + gc + lr*rsize + offset));
  vstore8(tmptot1, 0, &tot.t[lr * LOSI]);
  tmptot1 = convert_float8(vload8(0, pic + gc + 8 + lr*rsize + offset));
  vstore8(tmptot1, 0, &tot.t[lr * LOSI + DCTSIZE]);

  barrier(CLK_LOCAL_MEM_FENCE);

  jpeg_fdct_float(ctblock1, &tot.t[lr], LOSI);
#pragma unroll
  for (int i = 0; i < DCTSIZE2; i++) {
    ctblock1[i] = round(ctblock1[i] * fquant[i]) * iquant[i];
  }

  jpeg_idct_float(ctblock1, ctblock1, DCTSIZE);

  barrier(CLK_LOCAL_MEM_FENCE);
  vstore16((float16)(0), 0, &tot.f[lr * 16]);
  barrier(CLK_LOCAL_MEM_FENCE);

#pragma unroll
  for (int y = 0; y < DCTSIZE; y++)
  {
#pragma unroll
    for (int z = 0; z < DCTSIZE; z++)
    {
      tot.f[lr + z + y * LOSI] += ctblock1[y * 8 + z];
      barrier(CLK_LOCAL_MEM_FENCE);
    }
  }

  //ushort8 tmptotfin1 = vload8(0, res + gc + lr*rsize + offset) + convert_ushort8_sat_rtz(vload8(0, &tot.f[lr*LOSI]));
  //vstore8(tmptotfin1, 0, res + gc + lr*rsize + offset);
  //tmptotfin1 = vload8(0, res + gc + 8 + lr*rsize + offset) + convert_ushort8_sat_rtz(vload8(0, &tot.f[lr * LOSI + DCTSIZE]));
  //vstore8(tmptotfin1, 0, res + gc + 8 + lr*rsize + offset);
  *(__global ushort8 *)(res + gc + lr*rsize + offset) += convert_ushort8_sat_rtz(*(__local float8 *)&tot.f[lr*LOSI]);
  *(__global ushort8 *)(res + gc + 8 + lr*rsize + offset) += convert_ushort8_sat_rtz(*(__local float8 *)&tot.f[lr * LOSI + DCTSIZE]);

}

__kernel void devultur_nook(const __global uchar* pic, __global ushort* res, const __global float *fquant, const __global float *iquant, int rsize, int csize, int offset)
{
  uint glcrt = get_local_id(0);
  __local uchar input_buf[DCTSIZE2];
  __local float dataptr[DCTSIZE2];
  uchar8 tmptot = vload8(0, pic + glcrt*rsize + offset);
  vstore8(tmptot, 0, &input_buf[glcrt * 8]);
  barrier(CLK_LOCAL_MEM_FENCE);

  float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  float tmp10, tmp11, tmp12, tmp13;
  float z1, z2, z3, z4, z5, z11, z13;
  float z10, z12;

  /* pass 1: process rows. */
  //dataptr = outptr;
  tmp0 = convert_float(input_buf[glcrt * 8 + 0] + input_buf[glcrt * 8 + 7]);
  tmp7 = convert_float(input_buf[glcrt * 8 + 0] - input_buf[glcrt * 8 + 7]);
  tmp1 = convert_float(input_buf[glcrt * 8 + 1] + input_buf[glcrt * 8 + 6]);
  tmp6 = convert_float(input_buf[glcrt * 8 + 1] - input_buf[glcrt * 8 + 6]);
  tmp2 = convert_float(input_buf[glcrt * 8 + 2] + input_buf[glcrt * 8 + 5]);
  tmp5 = convert_float(input_buf[glcrt * 8 + 2] - input_buf[glcrt * 8 + 5]);
  tmp3 = convert_float(input_buf[glcrt * 8 + 3] + input_buf[glcrt * 8 + 4]);
  tmp4 = convert_float(input_buf[glcrt * 8 + 3] - input_buf[glcrt * 8 + 4]);

  tmp10 = tmp0 + tmp3;
  tmp13 = tmp0 - tmp3;
  tmp11 = tmp1 + tmp2;
  tmp12 = tmp1 - tmp2;

  /* apply unsigned->signed conversion. */
  dataptr[glcrt * 8 + 0] = tmp10 + tmp11 - 8 * CENTERJSAMPLE;
  dataptr[glcrt * 8 + 4] = tmp10 - tmp11;

  z1 = (tmp12 + tmp13) * 0.707106781f;
  dataptr[glcrt * 8 + 2] = tmp13 + z1;
  dataptr[glcrt * 8 + 6] = tmp13 - z1;

  tmp10 = tmp4 + tmp5;
  tmp11 = tmp5 + tmp6;
  tmp12 = tmp6 + tmp7;

  z5 = (tmp10 - tmp12) * 0.382683433f;
  z2 = 0.541196100f * tmp10 + z5;
  z4 = 1.306562965f * tmp12 + z5;
  z3 = tmp11 * 0.707106781f;

  z11 = tmp7 + z3;
  z13 = tmp7 - z3;

  dataptr[glcrt * 8 + 5] = z13 + z2;
  dataptr[glcrt * 8 + 3] = z13 - z2;
  dataptr[glcrt * 8 + 1] = z11 + z4;
  dataptr[glcrt * 8 + 7] = z11 - z4;

  //dataptr += DCTSIZE;
  //   input_buf += input_stride;
  barrier(CLK_LOCAL_MEM_FENCE);
  /* pass 2: process columns. */

  //dataptr = outptr;
  //for (ctr = DCTSIZE - 1; ctr >= 0; ctr--) {
  tmp0 = dataptr[glcrt + DCTSIZE * 0] + dataptr[glcrt + DCTSIZE * 7];
  tmp7 = dataptr[glcrt + DCTSIZE * 0] - dataptr[glcrt + DCTSIZE * 7];
  tmp1 = dataptr[glcrt + DCTSIZE * 1] + dataptr[glcrt + DCTSIZE * 6];
  tmp6 = dataptr[glcrt + DCTSIZE * 1] - dataptr[glcrt + DCTSIZE * 6];
  tmp2 = dataptr[glcrt + DCTSIZE * 2] + dataptr[glcrt + DCTSIZE * 5];
  tmp5 = dataptr[glcrt + DCTSIZE * 2] - dataptr[glcrt + DCTSIZE * 5];
  tmp3 = dataptr[glcrt + DCTSIZE * 3] + dataptr[glcrt + DCTSIZE * 4];
  tmp4 = dataptr[glcrt + DCTSIZE * 3] - dataptr[glcrt + DCTSIZE * 4];

  tmp10 = tmp0 + tmp3;
  tmp13 = tmp0 - tmp3;
  tmp11 = tmp1 + tmp2;
  tmp12 = tmp1 - tmp2;

  dataptr[glcrt + DCTSIZE * 0] = tmp10 + tmp11;
  dataptr[glcrt + DCTSIZE * 4] = tmp10 - tmp11;

  z1 = (tmp12 + tmp13) * 0.707106781f;
  dataptr[glcrt + DCTSIZE * 2] = tmp13 + z1;
  dataptr[glcrt + DCTSIZE * 6] = tmp13 - z1;

  tmp10 = tmp4 + tmp5;
  tmp11 = tmp5 + tmp6;
  tmp12 = tmp6 + tmp7;

  z5 = (tmp10 - tmp12) * 0.382683433f;
  z2 = 0.541196100f * tmp10 + z5;
  z4 = 1.306562965f * tmp12 + z5;
  z3 = tmp11 * 0.707106781f;

  z11 = tmp7 + z3;
  z13 = tmp7 - z3;

  dataptr[glcrt + DCTSIZE * 5] = z13 + z2;
  dataptr[glcrt + DCTSIZE * 3] = z13 - z2;
  dataptr[glcrt + DCTSIZE * 1] = z11 + z4;
  dataptr[glcrt + DCTSIZE * 7] = z11 - z4;

#pragma unroll
  for (int i = 0; i < DCTSIZE; i++) {
    dataptr[glcrt + i * 8] = round(dataptr[glcrt + i * 8] * fquant[glcrt + i * 8]) * iquant[glcrt + i * 8];
  }

  tmp0 = dataptr[glcrt + DCTSIZE * 0];
  tmp1 = dataptr[glcrt + DCTSIZE * 2];
  tmp2 = dataptr[glcrt + DCTSIZE * 4];
  tmp3 = dataptr[glcrt + DCTSIZE * 6];

  tmp10 = tmp0 + tmp2;
  tmp11 = tmp0 - tmp2;

  tmp13 = tmp1 + tmp3;
  tmp12 = (tmp1 - tmp3) * 1.414213562f - tmp13;

  tmp0 = tmp10 + tmp13;
  tmp3 = tmp10 - tmp13;
  tmp1 = tmp11 + tmp12;
  tmp2 = tmp11 - tmp12;

  tmp4 = dataptr[glcrt + DCTSIZE * 1];
  tmp5 = dataptr[glcrt + DCTSIZE * 3];
  tmp6 = dataptr[glcrt + DCTSIZE * 5];
  tmp7 = dataptr[glcrt + DCTSIZE * 7];

  z13 = tmp6 + tmp5;
  z10 = tmp6 - tmp5;
  z11 = tmp4 + tmp7;
  z12 = tmp4 - tmp7;

  tmp7 = z11 + z13;
  tmp11 = (z11 - z13) * 1.414213562f;

  z5 = (z10 + z12) * 1.847759065f;
  tmp10 = z5 - z12 * 1.082392200f;
  tmp12 = z5 - z10 * 2.613125930f;

  tmp6 = tmp12 - tmp7;
  tmp5 = tmp11 - tmp6;
  tmp4 = tmp10 - tmp5;

  dataptr[glcrt + DCTSIZE * 0] = tmp0 + tmp7;
  dataptr[glcrt + DCTSIZE * 7] = tmp0 - tmp7;
  dataptr[glcrt + DCTSIZE * 1] = tmp1 + tmp6;
  dataptr[glcrt + DCTSIZE * 6] = tmp1 - tmp6;
  dataptr[glcrt + DCTSIZE * 2] = tmp2 + tmp5;
  dataptr[glcrt + DCTSIZE * 5] = tmp2 - tmp5;
  dataptr[glcrt + DCTSIZE * 3] = tmp3 + tmp4;
  dataptr[glcrt + DCTSIZE * 4] = tmp3 - tmp4;

  /* pass 2: process rows. */
  barrier(CLK_LOCAL_MEM_FENCE);
  /* prepare range-limit and float->int conversion */
  z5 = dataptr[glcrt * 8 + 0] + (CENTERJSAMPLE + 0.5f);
  tmp10 = z5 + dataptr[glcrt * 8 + 4];
  tmp11 = z5 - dataptr[glcrt * 8 + 4];

  tmp13 = dataptr[glcrt * 8 + 2] + dataptr[glcrt * 8 + 6];
  tmp12 = (dataptr[glcrt * 8 + 2] - dataptr[glcrt * 8 + 6]) * 1.414213562f - tmp13;

  tmp0 = tmp10 + tmp13;
  tmp3 = tmp10 - tmp13;
  tmp1 = tmp11 + tmp12;
  tmp2 = tmp11 - tmp12;

  z13 = dataptr[glcrt * 8 + 5] + dataptr[glcrt * 8 + 3];
  z10 = dataptr[glcrt * 8 + 5] - dataptr[glcrt * 8 + 3];
  z11 = dataptr[glcrt * 8 + 1] + dataptr[glcrt * 8 + 7];
  z12 = dataptr[glcrt * 8 + 1] - dataptr[glcrt * 8 + 7];

  tmp7 = z11 + z13;
  tmp11 = (z11 - z13) * 1.414213562f;

  z5 = (z10 + z12) * 1.847759065f;
  tmp10 = z5 - z12 * 1.082392200f;
  tmp12 = z5 - z10 * 2.613125930f;

  tmp6 = tmp12 - tmp7;
  tmp5 = tmp11 - tmp6;
  tmp4 = tmp10 - tmp5;

  /* final output stage: float->int conversion and range-limit */
  input_buf[glcrt * 8 + 0] = range_limit(tmp0 + tmp7);
  input_buf[glcrt * 8 + 7] = range_limit(tmp0 - tmp7);
  input_buf[glcrt * 8 + 1] = range_limit(tmp1 + tmp6);
  input_buf[glcrt * 8 + 6] = range_limit(tmp1 - tmp6);
  input_buf[glcrt * 8 + 2] = range_limit(tmp2 + tmp5);
  input_buf[glcrt * 8 + 5] = range_limit(tmp2 - tmp5);
  input_buf[glcrt * 8 + 3] = range_limit(tmp3 + tmp4);
  input_buf[glcrt * 8 + 4] = range_limit(tmp3 - tmp4);
  barrier(CLK_LOCAL_MEM_FENCE);

  //ushort8 tmptotfin = vload8(0, res + glcrt*rsize + offset) + convert_ushort8_sat_rtz(vload8(0, &input_buf[glcrt * 8]));
  //vstore8(tmptotfin, 0, res + glcrt*rsize + offset);

  *(__global ushort8 *)(res + glcrt*rsize + offset) += convert_ushort8_sat_rtz(*(__local uchar8 *)&input_buf[glcrt * 8]);
}

__kernel void devultur_division(__global ushort* res, int w, int h)
{
  uint j = get_global_id(1);
  uint i = get_global_id(0);
  int ky = j < DCTSIZE ? j + 1 : j > h - DCTSIZE ? h - j : DCTSIZE;
  int kx = i < DCTSIZE ? i + 1 : i > w - DCTSIZE ? w - i : DCTSIZE;
  res[i + j*w] /= kx*ky;
}

__kernel void devultur_devnull(__global ushort2* res)
{
  res[get_global_id(0)] = 0;
}