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

    conf.setenv("cenv")
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
       if not conf.env.LIB_FFTW3:
          # Try to link the fftw without any further options.
          conf.check(lib='fftw3', uselib_store='FFTW3', mandatory=False)

    if conf.env.LIB_FFTW3:
       conf.all_envs[''].FCFLAGS_FFTW3 = conf.env.CFLAGS_FFTW3
       conf.all_envs[''].LIB_FFTW3 = conf.env.LIB_FFTW3
       conf.all_envs[''].LIBPATH_FFTW3 = conf.env.LIBPATH_FFTW3
       conf.all_envs[''].INCLUDES_FFTW3 = conf.env.INCLUDES_FFTW3
    conf.setenv('')

    if conf.env.LIB_FFTW3:
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

    # Check for m lib (necessary for FXTP):
    conf.check_cc(lib='m', uselib_store='MATH')



def build(bld):
    from waflib import Logs
    from waflib.extras.utest_results import utests

    ply_sources = bld.path.ant_glob('source/fpt/*.f90')
    ply_sources += bld.path.ant_glob('source/*.f90')
    ply_ppsources = bld.path.ant_glob('source/fpt/*.fpp')
    ply_ppsources += bld.path.ant_glob('source/*.fpp')
    ply_sources += ply_ppsources

    #FXTP Sources
    fxtp_wrap_sources = ['external/fxtp/fxt_fif.f90',
                         'external/fxtp/fxt_fwrap.f90']

    fxtp_sources = ['fxt_faltld.c',
		'fxt_faltld_comp.c',
		'fxt_faltld_preproc.c',
		'fxt_file.c',
		'fxt_error.c',
                'fxt_flptld.c',
		'fxt_flptld_comp.c', 
		'fxt_flptld_preproc.c',
		'fxt_fmmld.c',
		'fxt_fmmld_estimate.c',
		'fxt_fmmld_evaluate.c',
		'fxt_fmmld_evl.c',
		'fxt_fmmld_exp.c',
		'fxt_fmmld_lagld.c',
		'fxt_fmmld_locexpmat.c',
		'fxt_fmmld_matrix.c',
		'fxt_fmmld_mulexpmat.c',
		'fxt_fmmld_network.c',
		'fxt_fmmld_preproc.c',
		'fxt_fmmld_regions.c',
		'fxt_fmmld_scalefmm.c',
		'fxt_fmmld_scalevec.c',
		'fxt_fxtld.c',
		'fxt_fxtld_amat.c',
		'fxt_fxtld_comp.c',
		'fxt_fxtld_div.c',
		'fxt_fxtld_fmm.c',
		'fxt_fxtld_info.c',
		'fxt_fxtld_lagld.c',
		'fxt_fxtld_make.c',
		'fxt_fxtld_scale.c',
		'fxt_lagld.c',
		'fxt_lagld_comp.c',
		'fxt_lagld_drop.c',
		'fxt_lagld_normal.c',
		'fxt_lagld_oper.c',
		'fxt_math_gaussll.c',
		'fxt_math_legendrell.c',
		'fxt_math_lsqrt.c',
		'fxt_matld.c',
		'fxt_matld_drop.c',
		'fxt_matld_glq.c',
		'fxt_matld_ipstab.c',
		'fxt_matld_lagld.c',
		'fxt_matld_lowrank.c',
		'fxt_matld_vecld.c',
		'fxt_matll.c',
		'fxt_matll_drop.c',
		'fxt_matll_ipstab.c',
		'fxt_sarld.c',
		'fxt_sarld_comp.c',
		'fxt_sarld_error.c',
		'fxt_sarld_lagld.c',
		'fxt_sarld_scale.c',
		'fxt_vecl.c',
		'fxt_vecld.c',
		'fxt_vecll.c']
    
    for i_source in range(0, len(fxtp_sources)):
       fxtp_sources[i_source] = 'external/fxtp/fxtpack140715/' + fxtp_sources[i_source]


    if bld.cmd != 'gendoxy':
       if bld.env.LIB_FFTW3:
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
  
       bld( features = 'c',
            source = fxtp_sources,
            use = ['MATH'],
            includes = 'external/fxtp/fxtpack140715',
            target = 'fxtp_obj')

       bld( features = 'c',
            source = ['external/fxtp/fxtf_wrapper.c'],
            use = ['MATH'],
            includes = 'external/fxtp/fxtpack140715',
            target = 'fxtp_wrapper')

       bld( features = 'fc',
            source = fxtp_wrap_sources,
            use = ['MATH'],
            target = 'fxtp_wrap_obj')

       bld(
           features = 'coco fc',
           source = ply_sources,
           use = ['FFTW3', 'NAG', 'tem_objs', 'fftw_mod_obj', 'aotus'],
           target = 'ply_objs')

       bld(
           features = 'fc fcprogram',
           source = ['peons/approximate_1D_jump.f90'],
           use = ['FFTW3', 'NAG', 'tem_objs', 'ply_objs', 'fftw_mod_obj',
                  bld.env.mpi_mem_c_obj, 'fxtp_wrap_obj', 'fxtp_obj',
                  'PRECICE', 'fxtp_wrapper', 'aotus'],
           target = 'approximate_1D_jump')

       test_dep = ['FFTW3', 'NAG', 'tem_objs', bld.env.mpi_mem_c_obj,
                   'ply_objs', 'fftw_mod_obj', 'fxtp_wrap_obj', 'fxtp_obj',
                   'fxtp_wrapper', 'aotus']

       utest_sources = bld.path.ant_glob('utests/*_module.f90')

       bld(
         features = 'fc',
         source   = utest_sources,
         use      = ['tem_objs', 'aotus', 'ply_objs'],
         target   = 'ply_utest_objs')

       test_dep.append('ply_utest_objs')
       test_dep.append('PRECICE')

       utests(bld = bld, use = test_dep)

       if bld.env.LIB_FFTW3:
          utests(bld = bld, use = test_dep, path = 'utests/with_fftw')


    else:
       bld(
           features = 'coco',
           source   = ply_ppsources)
