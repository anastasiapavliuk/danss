#include "writer.h"

void write_x2(int ntimes) {
    for (int i = 0; i<ntimes; i++) {
        int k = 0;
        std::cin >> k;
        std::cout << k << " x 2 = " << k*2 << std::endl;
    }
}
