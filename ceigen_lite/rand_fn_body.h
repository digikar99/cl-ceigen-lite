

#define rand_fn_body(name, eigen_name, xtype, eigen_type, urng, ...) \
  extern "C"{                                                           \
    void CEIGEN_LITE_##name(const long n, xtype *x, __VA_ARGS__){       \
      Map<eigen_type> eigmat(x, n, 1);                                  \
      eigmat = eigen_name<eigen_type>(n, 1, urng, );                     \
      return;                                                           \
    }                                                                   \
  };
