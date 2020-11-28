
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
void f6print(uint64_t* p){
  f2print(p);
  f2print(p+12);
  f2print(p+24);
}
void f12print(uint64_t* p){
  f6print(p);
  f6print(p+36);
}


void test1(){
  // using src/e2.c and src/e1.c, negative generators [in Montgomery]
  blst_p1_affine BLS12_381_NEG_G1 = { /* negative generator [in Montgomery] */
  /* (0x17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905
   *    a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb << 384) % P */
  { TO_LIMB_T(0x5cb38790fd530c16), TO_LIMB_T(0x7817fc679976fff5),
    TO_LIMB_T(0x154f95c7143ba1c1), TO_LIMB_T(0xf0ae6acdf3d0e747),
    TO_LIMB_T(0xedce6ecc21dbf440), TO_LIMB_T(0x120177419e0bfb75) },
  /* (0x114d1d6855d545a8aa7d76c8cf2e21f267816aef1db507c9
   *    6655b9d5caac42364e6f38ba0ecb751bad54dcd6b939c2ca << 384) % P */
  { TO_LIMB_T(0xff526c2af318883a), TO_LIMB_T(0x92899ce4383b0270),
    TO_LIMB_T(0x89d7738d9fa9d055), TO_LIMB_T(0x12caf35ba344c12a),
    TO_LIMB_T(0x3cff1b76964b5317), TO_LIMB_T(0x0e44d2ede9774430) },
};

  blst_p2_affine BLS12_381_NEG_G2 = { /* negative generator [in Montgomery] */
{ /* (0x024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02
        b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8 << 384) % P */
  { TO_LIMB_T(0xf5f28fa202940a10), TO_LIMB_T(0xb3f5fb2687b4961a),
    TO_LIMB_T(0xa1a893b53e2ae580), TO_LIMB_T(0x9894999d1a3caee9),
    TO_LIMB_T(0x6f67b7631863366b), TO_LIMB_T(0x058191924350bcd7) ,
  /* (0x13e02b6052719f607dacd3a088274f65596bd0d09920b61a
        b5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e << 384) % P */
   TO_LIMB_T(0xa5a9c0759e23f606), TO_LIMB_T(0xaaa0c59dbccd60c3),
    TO_LIMB_T(0x3bb17e18e2867806), TO_LIMB_T(0x1b1ab6cc8541b367),
    TO_LIMB_T(0xc2b6ed0ef2158547), TO_LIMB_T(0x11922a097360edf3) }
},
{ /* (0x0d1b3cc2c7027888be51d9ef691d77bcb679afda66c73f17
        f9ee3837a55024f78c71363275a75d75d86bab79f74782aa << 384) % P */
  { TO_LIMB_T(0x6d8bf5079fb65e61), TO_LIMB_T(0xc52f05df531d63a5),
    TO_LIMB_T(0x7f4a4d344ca692c9), TO_LIMB_T(0xa887959b8577c95f),
    TO_LIMB_T(0x4347fe40525c8734), TO_LIMB_T(0x197d145bbaff0bb5) ,
  /* (0x13fa4d4a0ad8b1ce186ed5061789213d993923066dddaf10
        40bc3ff59f825c78df74f2d75467e25e0f55f8a00fa030ed << 384) % P */
   TO_LIMB_T(0x0c3e036d209afa4e), TO_LIMB_T(0x0601d8f4863f9e23),
    TO_LIMB_T(0xe0832636bacc0a84), TO_LIMB_T(0xeb2def362a476f84),
    TO_LIMB_T(0x64044f659f0ee1e9), TO_LIMB_T(0x0ed54f48d5a1caa7) }
}
};


  blst_fp12 *myfp12 = malloc(48*12);
  blst_p2_affine* myp2 = &BLS12_381_NEG_G2;
  blst_p1_affine* myp1 = &BLS12_381_NEG_G1;

  blst_miller_loop(myfp12,myp2,myp1);
}

void test2(){
  // identity elements, copied from https://github.com/ethereum/EIPs/blob/master/EIPS/eip-2539.md#specification
  uint8_t P[2*48];
  uint8_t Q[4*48];
  uint8_t* P0 = P;
  uint8_t* P1 = P+48;
  uint8_t* Q0 = Q;
  uint8_t* Q1 = Q0+48;
  uint8_t* Q2 = Q1+48;
  uint8_t* Q3 = Q2+48;
  hexstr_to_bytearray(P0,"008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef");
  hexstr_to_bytearray(P1,"01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6");
  hexstr_to_bytearray(Q0,"018480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196");
  hexstr_to_bytearray(Q1,"00ea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe");
  hexstr_to_bytearray(Q2,"00690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf");
  hexstr_to_bytearray(Q3,"00f8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93");

  blst_fp12 *myfp12 = malloc(48*12);

  blst_miller_loop(myfp12,(blst_p2_affine*)Q,(blst_p1_affine*)P);

}

void test2_full(){
  // identity elements, copied from https://github.com/ethereum/EIPs/blob/master/EIPS/eip-2539.md#specification
  uint8_t P[2*48];
  uint8_t Q[4*48];
  uint8_t* P0 = P;
  uint8_t* P1 = P+48;
  uint8_t* Q0 = Q;
  uint8_t* Q1 = Q0+48;
  uint8_t* Q2 = Q1+48;
  uint8_t* Q3 = Q2+48;
  hexstr_to_bytearray(P0,"008848defe740a67c8fc6225bf87ff5485951e2caa9d41bb188282c8bd37cb5cd5481512ffcd394eeab9b16eb21be9ef");
  hexstr_to_bytearray(P1,"01914a69c5102eff1f674f5d30afeec4bd7fb348ca3e52d96d182ad44fb82305c2fe3d3634a9591afd82de55559c8ea6");
  hexstr_to_bytearray(Q0,"018480be71c785fec89630a2a3841d01c565f071203e50317ea501f557db6b9b71889f52bb53540274e3e48f7c005196");
  hexstr_to_bytearray(Q1,"00ea6040e700403170dc5a51b1b140d5532777ee6651cecbe7223ece0799c9de5cf89984bff76fe6b26bfefa6ea16afe");
  hexstr_to_bytearray(Q2,"00690d665d446f7bd960736bcbb2efb4de03ed7274b49a58e458c282f832d204f2cf88886d8c7c2ef094094409fd4ddf");
  hexstr_to_bytearray(Q3,"00f8169fd28355189e549da3151a70aa61ef11ac3d591bf12463b01acee304c24279b83f5e52270bd9a1cdd185eb8f93");

  blst_fp12 *myfp12 = malloc(48*12);

  blst_miller_loop(myfp12,(blst_p2_affine*)Q,(blst_p1_affine*)P);
  blst_fp12 *myfp12_out = malloc(48*12);
  blst_final_exp(myfp12_out, myfp12);
  f12print((uint64_t*)myfp12_out);
}


bool pairing_eq2(uint8_t *pG1_1, uint8_t *pG2_1, uint8_t *pG1_2, uint8_t *pG2_2) {
/*
2pairing check:
blst_miller_loop()
blst_miller_loop()
blst_fp12_mul()
blst_final_exp()
blst_fp12_is_one()
*/
    uint8_t output_miller1[48 * 12];
    uint8_t output_miller2[48 * 12];
    uint8_t final_exp_input[48 * 12];
    uint8_t output[48 * 12];

    blst_miller_loop(output_miller1, pG2_1, pG1_1);
    blst_miller_loop(output_miller2, pG2_2, pG1_2);
    blst_fp12_mul(final_exp_input, output_miller1, output_miller2);
    blst_final_exp(result, final_exp_input);
    return blst_is_one(result);
}

// test e(pG1, pG2) * e(-pG1, pG2) == 1 where pG1 and pG2 are generator points
// reference: {{websnark branch link}}
void test_pairing2_check() {
uint8_t p1_G1[48], p1_G2[48 * 2], p2_G1[48], p2_n_G1[48 * 2];

// p1_G1 <- G1 generator (normal form)
// p1_G2 <- G2 generator (normal form)

// p2_G1 <- neg(p1_G1)
// p2_G2 <- p1_G2


/*
assert e(p1[0], p2[0]) * e(p1[1], p2[1]) == 1
*/
}

void test3(){	// from wasmsnark/test
// uncommented miller loop inputs/output in test/bls12381.js
//cd wasmsnark && ~/repos/node/node-v12.18.4-linux-x64/bin/npx mocha test/bls12381.js

  uint8_t P[2*48];
  uint8_t Q[4*48];
  uint8_t* P0 = P;
  uint8_t* P1 = P+48;
  uint8_t* Q0 = Q;
  uint8_t* Q1 = Q0+48;
  uint8_t* Q2 = Q1+48;
  uint8_t* Q3 = Q2+48;
  hexstr_to_bytearray(P0,"0f81da25ecf1c84b577fefbedd61077a81dc43b00304015b2b596ab67f00e41c86bb00ebd0f90d4b125eb0539891aeed");
  hexstr_to_bytearray(P1,"11af629591ec86916d6ce37877b743fe209a3af61147996c1df7fd1c47b03181cd806fd31c3071b739e4deb234bd9e19");
  hexstr_to_bytearray(Q0,"024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8");
  hexstr_to_bytearray(Q1,"13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e");
  hexstr_to_bytearray(Q2,"0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801");
  hexstr_to_bytearray(Q3,"0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be");

  blst_fp12 *myfp12 = malloc(48*12);

  blst_miller_loop(myfp12,(blst_p2_affine*)Q,(blst_p1_affine*)P);
  //f12print((uint64_t*)myfp12);

  //blst_fp12 *myfp12_out = malloc(48*12);
  //blst_final_exp(myfp12_out, myfp12);
  //f12print((uint64_t*)myfp12_out);

}

void test4(){
  // test from https://tools.ietf.org/id/draft-yonezawa-pairing-friendly-curves-02.html#rfc.appendix.B
  uint8_t P[2*48];
  uint8_t Q[4*48];
  uint8_t* P0 = P;
  uint8_t* P1 = P+48;
  uint8_t* Q0 = Q;
  uint8_t* Q1 = Q0+48;
  uint8_t* Q2 = Q1+48;
  uint8_t* Q3 = Q2+48;
  hexstr_to_bytearray(P0,"17f1d3a73197d7942695638c4fa9ac0fc3688c4f9774b905a14e3a3f171bac586c55e83ff97a1aeffb3af00adb22c6bb");
  hexstr_to_bytearray(P1,"08b3f481e3aaa0f1a09e30ed741d8ae4fcf5e095d5d00af600db18cb2c04b3edd03cc744a2888ae40caa232946c5e7e1");
  hexstr_to_bytearray(Q0,"024aa2b2f08f0a91260805272dc51051c6e47ad4fa403b02b4510b647ae3d1770bac0326a805bbefd48056c8c121bdb8");
  hexstr_to_bytearray(Q1,"13e02b6052719f607dacd3a088274f65596bd0d09920b61ab5da61bbdc7f5049334cf11213945d57e5ac7d055d042b7e");
  hexstr_to_bytearray(Q2,"0ce5d527727d6e118cc9cdc6da2e351aadfd9baa8cbdd3a76d429a695160d12c923ac9cc3baca289e193548608b82801");
  hexstr_to_bytearray(Q3,"0606c4a02ea734cc32acd2b02bc28b99cb3e287e85a763af267492ab572e99ab3f370d275cec1da1aaa9075ff05f79be");

  blst_fp12 *myfp12 = malloc(48*12);

  blst_miller_loop(myfp12,(blst_p2_affine*)Q,(blst_p1_affine*)P);
  //printf("output of miller loop\n");
  //f12print((uint64_t*)myfp12);

  blst_fp12 *myfp12_out = malloc(48*12);
  blst_final_exp(myfp12_out, myfp12);
  f12print((uint64_t*)myfp12_out);
/* output in spec is:
      0x11619b45f61edfe3b47a15fac19442526ff489dcda25e59121d9931438907dfd448299a87dde3a649bdba96e84d54558
      0x153ce14a76a53e205ba8f275ef1137c56a566f638b52d34ba3bf3bf22f277d70f76316218c0dfd583a394b8448d2be7f
      0x095668fb4a02fe930ed44767834c915b283b1c6ca98c047bd4c272e9ac3f3ba6ff0b05a93e59c71fba77bce995f04692
      0x16deedaa683124fe7260085184d88f7d036b86f53bb5b7f1fc5e248814782065413e7d958d17960109ea006b2afdeb5f
      0x09c92cf02f3cd3d2f9d34bc44eee0dd50314ed44ca5d30ce6a9ec0539be7a86b121edc61839ccc908c4bdde256cd6048
      0x111061f398efc2a97ff825b04d21089e24fd8b93a47e41e60eae7e9b2a38d54fa4dedced0811c34ce528781ab9e929c7
      0x01ecfcf31c86257ab00b4709c33f1c9c4e007659dd5ffc4a735192167ce197058cfb4c94225e7f1b6c26ad9ba68f63bc
      0x08890726743a1f94a8193a166800b7787744a8ad8e2f9365db76863e894b7a11d83f90d873567e9d645ccf725b32d26f
      0x0e61c752414ca5dfd258e9606bac08daec29b3e2c57062669556954fb227d3f1260eedf25446a086b0844bcd43646c10
      0x0fe63f185f56dd29150fc498bbeea78969e7e783043620db33f75a05a0a2ce5c442beaff9da195ff15164c00ab66bdde
      0x10900338a92ed0b47af211636f7cfdec717b7ee43900eee9b5fc24f0000c5874d4801372db478987691c566a8c474978
      0x1454814f3085f0e6602247671bc408bbce2007201536818c901dbd4d2095dd86c1ec8b888e59611f60a301af7776be3d
*/

}


int main(int argc,char**argv){



  //test1();
  //test2();
  test2_full();
  // test3();
  //test4();

  return 0;

}


