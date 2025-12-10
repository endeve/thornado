module ExternalTimersModule
   implicit none

   private

   abstract interface
      subroutine ExternalTimersType(name)
         implicit none
         character(len=*), intent(in) :: name
      end subroutine ExternalTimersType
   end interface

   procedure(ExternalTimersType), pointer, save :: ExternalTimerStart => ExternalTimer_Default
   procedure(ExternalTimersType), pointer, save :: ExternalTimerStop => ExternalTimer_Default

   public :: ExternalTimersType
   public :: SetExternalTimer
   public :: ExternalTimerStart
   public :: ExternalTimerStop

contains

   subroutine ExternalTimer_Default(name)
      character(len=*), intent(in) :: name
   end subroutine ExternalTimer_Default

   subroutine SetExternalTimer(TimerStart, TimerStop)
      procedure(ExternalTimersType) :: TimerStart, TimerStop

      ExternalTimerStart => TimerStart
      ExternalTimerStop  => TimerStop
   end subroutine SetExternalTimer
end module ExternalTimersModule
