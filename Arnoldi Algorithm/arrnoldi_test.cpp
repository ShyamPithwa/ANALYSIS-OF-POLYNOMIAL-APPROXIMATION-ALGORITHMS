#include "Arnoldi_iteration.h"
using namespace std;

int main(){
    Matrix data("dataset.txt");
    Arnoldi_Krylov_Polynomial polynomial(data);
    return 0;
}
