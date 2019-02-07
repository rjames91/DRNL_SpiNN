/*
 * IHC_AN.h
 *
 *  Created on: 1 Nov 2016
 *      Author: rjames
 */

#ifndef IHC_AN_softfloat_H_
#define IHC_AN_softfloat_H_

#define REAL double //float
#define REAL_CONST(x) x//x##f


#define MAX_CHIPX 1//255
#define MAX_CHIPY 1//255
#define MAX_COREID 16
//#define SEED_SEL_SIZE (((MAX_CHIPX << 8) | MAX_CHIPY) << 8) | (MAX_COREID << 3)   
#define SEED_SEL_SIZE 1024

#define SEGSIZE 8//16//96//100

#define SAMPLING_FREQUENCY 44100.
#define MOC_DELAY_MS 50

#define MAX_SIGNAL_S 1

#define TIMER2_CONF        0x82
#define TIMER2_LOAD        0

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define ABS(x) (((x)<0) ? -(x) : (x))
#define SIGN(x) (((x)<0) ? -(1.0) : (1.0))

typedef union
{
	uint32_t u;
	float f;
} uint_float_union;

typedef struct key_mask_table {
    uint32_t key;
    uint32_t mask;
    uint32_t conn_index;
} key_mask_table_entry;

typedef struct{
    uint32_t e_index;
    uint32_t w_index;
    uint32_t id_shift;
}last_neuron_info_t;

#endif /* IHC_AN_H_ */
