!< author: Arthur Francisco
!  version: 1.0.0
!  date: october, 23 2024
!
!  <span style="color: #337ab7; font-family: cabin; font-size: 1.2em;">
!        **Main program**
!  </span>

program main
!$ use omp_lib
use script, only : read_job
use crest_param
use data_arch, only : I4, R8

implicit none

   call prg_surf

   contains

   subroutine prg_surf()
   !================================================================================================
   !<@note Main function...
   !
   ! - retrieve script (job) file
   ! - read script
   ! - run specific functions associated to a script keyword
   !
   !<@endnote
   !------------------------------------------------------------------------------------------------
   implicit none

      character(len=128) :: arg_prg
      character(len=512) :: job_file
      character(len=008) :: chara_d
      character(len=010) :: chara_t
      integer(kind=I4)   :: var_i

      ! String initialisation
      arg_prg  = repeat(' ',len(arg_prg))
      JOB_FILE = repeat(' ',len(JOB_FILE))

      var_i = 1
      call get_command_argument(var_i, arg_prg)       ! argument one: argument string
      if (len_trim(arg_prg) == 0) then                ! if there is no job file, stop
         write(TER,*) 'no job file, stop'
         stop
      else
         job_file = trim(arg_prg)
      endif

      call read_job(job_file)   ! the program executes 'prg_repeat' times

      write(TER,*) 'Program completed'

   return
   endsubroutine prg_surf

endprogram main
