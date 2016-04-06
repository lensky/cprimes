#include <math.h>

#include "cvector.h"

void new_cvector_cap(cvector* cv, size_t sizeof_elt, size_t cap) {
    cv->length = 0;
    cv->capacity = cap;
    cv->sizeof_elt = sizeof_elt;

    cv->data = (cvector_data_unit*) malloc(sizeof_elt * cap);
}

void new_cvector(cvector* cv, size_t sizeof_elt) {
    return new_cvector_cap(cv, sizeof_elt, CVECTOR_DEF_CAP);
}

void free_cvector(cvector *cv) {
    free_cvector_data(cv->data);
}

void free_cvector_data(cvector_data_unit* data) {
    free(data);
}

cvector_data_unit* next_elt(cvector *cv) {
    if (cv->length == cv->capacity) {
        cv->capacity = (size_t) ceil((double) cv->capacity * CVECTOR_GROWTHF);
        cv->data = (cvector_data_unit*)
            realloc(cv->data, cv->capacity * cv->sizeof_elt);
    }
    // Note we increment cv->length below
    return cv->data + (cv->sizeof_elt * (cv->length++));
}
