//#include "mainwindow.h"
//#include <gtkmm/application.h>
#include<iostream>
#include <vector>       
#include "FS_make_kgrid.cpp"

int main(){
double SYMM_OP[18]={1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,1,1,1};
FS_make_kgrid(8,2, SYMM_OP);
}

