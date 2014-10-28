#! /usr/bin/env python
APPNAME = 'polynomials'
VERSION = '1'

top = '.'
out = 'build'

def configure(conf):
    from waflib import Errors, Logs

    # Check for FFTW lib
    fftw_libpath = None
    fftw_omplibpath = None

    if conf.options.fftw_libpath:
      fftw_libpath = conf.options.fftw_libpath
    elif conf.options.fftw_path:
      fftw_libpath = conf.options.fftw_path + '/lib'

    if fftw_libpath:
       Logs.info('Checking path to FFTW lib: ' + fftw_libpath)
       conf.check(lib='fftw3', libpath=fftw_libpath, uselib_store='FFTW3', mandatory=False)
       if conf.options.fftw_incpath:
         conf.env.INCLUDES_FFTW3 = conf.options.fftw_incpath
       elif conf.options.fftw_path:
         conf.env.INCLUDES_FFTW3 = conf.options.fftw_path+'/include'
    else:
       # Try to use pkg-config to find the FFTW library.
       conf.check_cfg(package='fftw3', uselib_store='FFTW3',
                      args=['--cflags', '--libs'], mandatory=False)
       if conf.env.LIB_FFTW3:
          conf.env.FCFLAGS_FFTW3 = conf.env.CFLAGS_FFTW3
       else:
          # Try to link the fftw without any further options.
          conf.check(lib='fftw3', uselib_store='FFTW3', mandatory=False)

    if conf.env.LIB_FFTW3:
       # Check for the OpenMP library of FFTW
       if conf.options.openmp:
         if conf.options.fftw_omplibpath:
           fftw_omplibpath = conf.options.fftw_omplibpath
         elif fftw_libpath:
           fftw_omplibpath = fftw_libpath

         if fftw_omplibpath:
           conf.check(lib='fftw3_omp', use='FFTW3', libpath=fftw_omplibpath, uselib_store='FFTW3_OMP', mandatory=True)
         else:
           conf.check(lib='fftw3_omp', use='FFTW3', uselib_store='FFTW3_OMP', mandatory=True)

       # Check for the fftw3.f03 header:
       try:
          conf.check_fc(fragment= "program test\n use, intrinsic :: iso_c_binding\n include 'fftw3.f03'\nend program test",
                        includes= conf.env.INCLUDES_FFTW3,
                        msg     = "Checking for the fftw3.f03 header",
                        errmsg  = "FFTW lib seems to be available, but no compatible Fortran header (need FFTW >=3.3.1)!"
                       )
       except Errors.ConfigurationError:
          conf.env.LIB_FFTW3 = None

    if not conf.env.LIB_FFTW3:
       # There is no FFTW found, use a dummy library to abort upon calls to the library.
       Logs.warn('There was NO FFTW3 library found, using a dummy implementation.')
       Logs.warn('You will not be able to use the FPT with the MODG scheme!')
       Logs.warn('')



def build(bld):
    from waflib import Logs

    ply_sources = bld.path.ant_glob('source/fpt/*.f90')
    ply_sources += bld.path.ant_glob('source/*.f90')
    ply_ppsources = bld.path.ant_glob('source/fpt/*.fpp')
    ply_ppsources += bld.path.ant_glob('source/*.fpp')
    ply_sources += ply_ppsources

    if bld.env.LIB_FFTW3:
       if bld.env.LIB_FFTW3_OMP:
          fftwdep = 'FFTW3_OMP FFTW3'
       else:
          fftwdep = 'FFTW3'
       bld( features = 'fc',
            source = 'external/fftw/fftw_wrap.f90',
            use = fftwdep,
            target = 'fftw_mod_obj')
    else:
       Logs.warn('Using the *dummy* FFTW here!')
       bld( features = 'fc',
            source = 'external/dummy/fftw_wrap.f90',
            target = 'fftw_mod_obj')

    bld(
        features = 'coco fc',
        source = ply_sources,
        use = ['FFTW3', 'NAG', 'tem_objs', 'fftw_mod_obj', 'aotus'],
        target = 'ply_objs')
