#include <stdio.h>
#include <stdlib.h>

void en_pot(double *posx, double *posy, double *posz, long ncharges, double *res) {
    double P = 0.0;
        for (long i = 0; i < ncharges; i++) {
            for (long j = 0; j < ncharges; j++) {
                if (i != j) {
                    double dx = posx[i] - posx[j];
                    double dy = posy[i] - posy[j];
                    double dz = posz[i] - posz[j];
                    double r = sqrt(dx*dx + dy*dy + dz*dz);
                    P += 1/r;
                }
            }
        }
        P *= 0.5;
        *res = P;
    }

