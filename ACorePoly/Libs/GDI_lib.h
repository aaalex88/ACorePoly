#ifndef GDI_LIB_H
#define GDI_LIB_H

#define WIN32_LEAN_AND_MEAN
#include <windows.h> 
#include <windowsx.h> 
#include <mmsystem.h> 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// DEFINES ////////////////////////////////////////////////
#define START 150
#define SCALE 0.01

// defines for windows 
#define WINDOW_CLASS_NAME "WINCLASS1"
#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 1000

// MACROS /////////////////////////////////////////////////

#define KEYDOWN(vk_code) ((GetAsyncKeyState(vk_code) & 0x8000) ? 1 : 0)
#define KEYUP(vk_code)   ((GetAsyncKeyState(vk_code) & 0x8000) ? 0 : 1)

// GLOBALS ////////////////////////////////////////////////
HWND      main_window_handle = NULL; // globally track main window

HINSTANCE hinstance_app      = NULL; // globally track hinstance
char buffer[80];                     // general printing buffer
WNDCLASSEX winclass; // this will hold the class we create
HWND	   hwnd;	 // generic window handle
MSG		   msg;		 // generic message
HDC        hdc;      // graphics device context



int gdi_index=0;





LRESULT CALLBACK WindowProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam);
int GDI_Init(HINSTANCE hInstance);
void GDI_Clear();
void GDI_DrawArray(int N, double * X, COLORREF col,  double scale=SCALE, int start=START);
void GDI_DrawArrayMod(int N, double * X, double * Y, COLORREF col,  double scale=SCALE, int start=START);
void GDI_DrawFormant(int i, double x, double scale=SCALE, double start=START);
void GDI_DrawACore(double fr, double step, int N, double * X, COLORREF col,  double scale=SCALE, int start=START);
//void GDI_Text(char* text, double x, double y);



#endif