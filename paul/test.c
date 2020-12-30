
/*
bash ../build.sh 
gcc test.c libblst.a -static
./a.out
*/

#include <stdio.h>
#include <stdlib.h>	// for malloc()
#include <string.h>	// for strlen()
#include "../bindings/blst.h"

#define TO_LIMB_T(limb64)	limb64

// hex string to int array conversion
// input is string of hex characters, without 0x prefix
// also converts to little endian (ie least significant nibble first)
void hexstr_to_bytearray(uint8_t* out, char* in){
  //printf("hexstr_to_intarray(%s)\n",in);
  size_t len = strlen(in);
  uint8_t byte = 0;
  uint8_t nibble = 0;
  int i;
  for (i=len-1; i>=0 ;i--){
    nibble = in[i];
    if (nibble >= '0' && nibble <= '9') nibble = nibble - '0';
    else if (nibble >= 'a' && nibble <='f') nibble = nibble - 'a' + 10;
    else if (nibble >= 'A' && nibble <='F') nibble = nibble - 'A' + 10;
    else printf("ERROR: %s is not hex string.\n",in);
    if (!((i+len%2)%2)) {
      byte = (nibble<<4) + byte;
      *(out+(len-i)/2-1) = byte;
      byte=0;
    }
    else byte = nibble;
  }
  if (byte)
    *(out+(len-i)/2-1) = byte;
}



void f1print(uint64_t* p){
  for (int i=0;i<6;i++)
    printf("%016lx ",(p)[i]);
  printf("\n");
}
void f2print(uint64_t* p){
  f1print(p);
  f1print(p+6);
}

void g2print(uint64_t *elem) {
  f2print(elem);
  f2print(elem + 12);
}

void f6print(uint64_t* p){
  f2print(p);
  f2print(p+12);
  f2print(p+24);
}
void f12print(uint64_t* p){
  f6print(p);
  f6print(p+36);
}

// test e(pG1_1, pG2_1) * e(pG1_2, pG2_2) == 1
bool pairing_eq2(blst_p1_affine *pG1_1, blst_p2_affine *pG2_1, blst_p1_affine *pG1_2, blst_p2_affine *pG2_2) {
  // not sure if blst cares about clobbering outputs.  non-clobbered here just to be safe
  uint8_t output_miller1[48 * 12];
  uint8_t output_miller2[48 * 12];
  uint8_t final_exp_input[48 * 12];
  uint8_t result[48 * 12];
  memcpy(final_exp_input, blst_fp12_one(), sizeof(blst_fp12));

    if(!blst_p1_affine_in_g1(pG1_1)) {
      printf("bad pG1_1\n");
      return 0;
    }

    if(!blst_p2_affine_in_g2(pG2_1)) {
      printf("bad pG2_1\n");
      return 0;
    }

    if(!blst_p1_affine_in_g1(pG1_2)) {
      printf("bad pG1_2\n");
      return 0;
    }

    if(!blst_p2_affine_in_g2(pG2_2)) {
      printf("bad pG2_2\n");
      return 0;
    }

  printf("pG1\n");
  f2print(pG1_1);

  printf("pG2\n");
  g2print(pG2_1);

  blst_miller_loop(output_miller1, pG2_1, pG1_1);
  blst_fp12_mul(final_exp_input, final_exp_input, output_miller1);

  blst_miller_loop(output_miller2, pG2_2, pG1_2);
  blst_fp12_mul(final_exp_input, final_exp_input, output_miller2);

  printf("final exp input\n");
  f12print(final_exp_input);

  blst_final_exp(result, final_exp_input);

  printf("final exp output\n");
  f12print(result);

  return blst_fp12_is_one(result);
}

bool pairing_eq1(blst_p1_affine *pG1_1, blst_p2_affine *pG2_1) {
    uint8_t input_final_exp[48 * 12];
    uint8_t result[48 * 12];

    if(!blst_p1_affine_in_g1(pG1_1)) {
      printf("bad pG1_1\n");
      return 0;
    }

    if(!blst_p2_affine_in_g2(pG2_1)) {
      printf("bad pG2_1\n");
      return 0;
    }

    blst_miller_loop(input_final_exp, pG2_1, pG1_1);
    //blst_fp12_mul(input_final_exp, output_miller1, output_miller2);

    blst_final_exp(result, input_final_exp);
    return blst_fp12_is_one(result);
}

// test e(pG1, pG2) * e(-pG1, pG2) == 1 where pG1 and pG2 are generator points
void test_pairing2_check_naive() {
  if (!pairing_eq2(&BLS12_381_G1, &BLS12_381_G2, &BLS12_381_NEG_G1, &BLS12_381_G2)) {
      printf("failed\n");
  }
}

// test e(pG1, pG2*42) * e(pG1, -pG2*42) == 1 where pG1 and pG2 are generator points
void test_pairing2_check() {
  blst_p1 pG1_1, pG1_2;
  blst_p2 pG2_1, pG2_2;

  printf("sitll runnling\n");

  blst_scalar s;
  uint32_t scalar_vals[8] = {0,0,0,0,0,0,42};
  blst_scalar_from_uint32(&s, scalar_vals);

  blst_p2_mult(&pG2_1, &BLS12_381_NEG_G2, &s, 192);
  blst_p2_mult(&pG2_2, &BLS12_381_G2, &s, 192);

  if (!pairing_eq2(&BLS12_381_G1, &pG2_1, &BLS12_381_G1, &pG2_2)) {
      printf("failed\n");
  }
}

int main(int argc,char**argv){
  test_pairing2_check();

  return 0;

}


