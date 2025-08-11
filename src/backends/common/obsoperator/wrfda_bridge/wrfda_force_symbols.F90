module wrfda_force_symbols
  use iso_c_binding
  use da_control, only: rootproc, kte
contains
  subroutine wrfda_force_ref() bind(C, name="wrfda_force_ref")
    integer :: i
    logical :: lr
    lr = rootproc
    i = kte
    if (lr) i = i + 1
  end subroutine wrfda_force_ref
end module wrfda_force_symbols


