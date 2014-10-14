#ifndef GDI_LIB_CPP
#define GDI_LIB_CPP



#include "GDI_lib.h"
#include "Acore_lib.h"




LRESULT CALLBACK WindowProc(HWND hwnd, UINT msg, WPARAM wparam, LPARAM lparam)
{
PAINTSTRUCT		ps;	
HDC				hdc;


switch(msg)
	{	
	case WM_CREATE: 
        {
		return(0);
		} break;
   
	case WM_PAINT: 
		{
		return(0);
   		} break;

	case WM_DESTROY: 
		{
		PostQuitMessage(0);
		return(0);
		} break;
	case WM_KEYDOWN:
		{
		} break;

	default:break;
    }
return (DefWindowProc(hwnd, msg, wparam, lparam));
}




int GDI_Init(HINSTANCE hInstance)
{
	winclass.cbSize         = sizeof(WNDCLASSEX);
	winclass.style			= CS_DBLCLKS | CS_OWNDC | 
		                      CS_HREDRAW | CS_VREDRAW;
	winclass.lpfnWndProc	= WindowProc;
	winclass.cbClsExtra		= 0;
	winclass.cbWndExtra		= 0;
	winclass.hInstance		= hInstance;
	winclass.hIcon			= LoadIcon(NULL, IDI_APPLICATION);
	winclass.hCursor		= LoadCursor(NULL, IDC_ARROW); 
	winclass.hbrBackground	= (HBRUSH)GetStockObject(BLACK_BRUSH);
	winclass.lpszMenuName	= NULL;
	winclass.lpszClassName	= (LPCWSTR)WINDOW_CLASS_NAME;
	winclass.hIconSm        = LoadIcon(NULL, IDI_APPLICATION);

	hinstance_app = hInstance;
	if (!RegisterClassEx(&winclass))
		return(0);
	if (!(hwnd = CreateWindowEx(NULL, (LPCWSTR)WINDOW_CLASS_NAME, (LPCWSTR)"DFT & Formants Demo", WS_OVERLAPPEDWINDOW | WS_VISIBLE,
					 	    0,0, WINDOW_WIDTH, WINDOW_HEIGHT, NULL, NULL, hInstance, NULL)))
	return(0);

	main_window_handle = hwnd;
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,RGB(200,200,200));
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
	MoveToEx(hdc, 0, START, NULL);
		LineTo(hdc,WINDOW_WIDTH, START);
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);
	return(1);
}



void GDI_Clear()
{
	HPEN hpen = CreatePen(PS_SOLID,1,RGB(0,0,0));
	HBRUSH hbr = CreateSolidBrush(RGB(0,0,0));
	SelectObject(hdc,hpen);
	SelectObject(hdc,hbr);
	Rectangle(hdc,0,0,WINDOW_WIDTH,WINDOW_HEIGHT);
}


void GDI_DrawArray(int N, double * X, COLORREF col,  double scale, int start)
{
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,col);
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
	MoveToEx(hdc, 0, start-scale*X[0], NULL);
	for (int i=1; i < N; i++)
	{
		LineTo(hdc,i, start-scale*X[i]);
	}
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);

	return;
}

void GDI_DrawArrayMod(int N, double * X, double * Y, COLORREF col,  double scale, int start)
{
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,col);
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
	MoveToEx(hdc, 0, start-scale*sqrt(X[0]*X[0]+Y[0]*Y[0]), NULL);
	for (int i=1; i < N; i++)
	{
		LineTo(hdc,i, start-scale*sqrt(X[i]*X[i]+Y[i]*Y[i]));
	}
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);

	return;
}

void GDI_DrawArrayEX(int N, double * X, COLORREF col,  double scale, int start)
{
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,col);
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
//	MoveToEx(hdc, 0, start-scale*X[0], NULL);
	for (int i=0; i < N; i++) // i=1
	{
		MoveToEx(hdc, i, start, NULL);
		LineTo(hdc,i, start-scale*X[i]);
//		SetPixel(hdc,i,start-scale*X[i],col);
//		LineTo(hdc,i, start-scale*X[i]);
	}
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);

	return;
}

void GDI_DrawArrayModEX(int N, double * X, double * Y, COLORREF col,  double scale, int start)
{
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,col);
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
//	MoveToEx(hdc, 0, start-scale*sqrt(X[0]*X[0]+Y[0]*Y[0]), NULL);
	for (int i=0; i < N; i++) // i=1
	{
		MoveToEx(hdc, i, start, NULL);
		LineTo(hdc,i, start-scale*sqrt(X[i]*X[i]+Y[i]*Y[i]));
//		SetPixel(hdc,i,start-scale*sqrt(X[i]*X[i]+Y[i]*Y[i]),col);
//		LineTo(hdc,i, start-scale*sqrt(X[i]*X[i]+Y[i]*Y[i]));
	}
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);

	return;
}

void GDI_DrawFormant(int i, double x, double scale, double start)
{
	hdc = GetDC(hwnd);
	SelectObject(hdc, (HGDIOBJ)WHITE_PEN);
    SelectObject(hdc, (HGDIOBJ)WHITE_BRUSH);
	Ellipse(hdc,i-5,start-scale*x-5,i+5,start-scale*x+5);
	ReleaseDC(hwnd, hdc);
}

void GDI_DrawACore(double fr, double step, int N, double * X, COLORREF col,  double scale, int start)
{
	if(fr<1.0) // fr==0.0
		return;
	hdc = GetDC(hwnd);
	HPEN hpen = CreatePen(PS_SOLID,1,col);
	HPEN old_hpen = (HPEN)SelectObject(hdc,hpen);
	
	for (int i=1; i*fr/step < N; i++)
	{
		MoveToEx(hdc, int(i*fr/step), start, NULL );
		LineTo(hdc,   int(i*fr/step), start-scale*X[int(i*fr/step)] );
	}
	SelectObject(hdc,old_hpen);
	DeleteObject(hpen);
	ReleaseDC(hwnd, hdc);

	return;
}

/*
void GDI_Text(char *text, double x, double y)
{
	if(!text)
		return;

	LPCTSTR t = TEXT(text);

	TextOut(hdc,x,y,t,strlen(text));
}
//*/



#endif