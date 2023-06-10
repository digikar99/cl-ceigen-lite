void CEIGEN_LITE_dlu(const long m, const long n,
                       double* a, char a_layout,
                       double* l, char l_layout,
                       double* u, char u_layout){
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, m, n);
      auto lu = amat.fullPivLu();
      if (l_layout == 'C' || l_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> lmat(l, m, m);
        lmat = lu.matrixLU().triangularView<Lower>();
        if (u_layout == 'C' || u_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }else if (u_layout == 'R' || u_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }
      }else if (l_layout == 'R' || l_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> lmat(l, m, m);
        lmat = lu.matrixLU().triangularView<Lower>();
        if (u_layout == 'C' || u_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }else if (u_layout == 'R' || u_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, m, n);
      auto lu = amat.fullPivLu();
      if (l_layout == 'C' || l_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> lmat(l, m, m);
        lmat = lu.matrixLU().triangularView<Lower>();
        if (u_layout == 'C' || u_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }else if (u_layout == 'R' || u_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }
      }else if (l_layout == 'R' || l_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> lmat(l, m, m);
        lmat = lu.matrixLU().triangularView<Lower>();
        if (u_layout == 'C' || u_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }else if (u_layout == 'R' || u_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> umat(u, m, n);
          umat = lu.matrixLU().triangularView<Upper>();
        }
      }
    }
  }
