#ifndef CVECTOR_H
#define CVECTOR_H

#include <stdint.h>
#include <stdlib.h>

#define CVECTOR_DEF_CAP 1024
#define CVECTOR_GROWTHF 1.5

typedef char cvector_data_unit; // sizeof(cvector_data_unit) == 1

typedef struct cvector {
    cvector_data_unit* data;
    size_t length;
    size_t capacity;
    size_t sizeof_elt;
} cvector;

void new_cvector_cap(cvector* cv, size_t sizeof_elt, size_t cap);
void new_cvector(cvector* cv, size_t sizeof_elt);
void free_cvector(cvector *cv);
void free_cvector_data(cvector_data_unit* data);

cvector_data_unit* next_elt(cvector *cv);

#endif /* CVECTOR_H */
