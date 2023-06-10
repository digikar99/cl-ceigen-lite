

#define one_arg_fn_body(name, eigen_name, itype, eigen_itype, otype, eigen_otype) \
  extern "C"{                                                           \
  void CEIGEN_LITE_##name(const long n, itype *x, otype *y){        \
    Map<eigen_itype> in(x, n, 1);                                   \
    Map<eigen_otype> out(y, n, 1);                                  \
    out = in.eigen_name();                                          \
    return;                                                         \
  }                                                                 \
  };
