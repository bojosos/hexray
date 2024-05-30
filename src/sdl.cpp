/***************************************************************************
 *   Copyright (C) 2009-2024 by Veselin Georgiev, Slavomir Kaslev,         *
 *                              Deyan Hadzhiev et al                       *
 *   admin@raytracing-bg.net                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/**
 * @File sdl.cpp
 * @Brief Implements the interface to SDL (mainly drawing to screen functions)
 */
#include <SDL.h>
#include <SDL_video.h>
#include <stdio.h>
#include "sdl.h"
#include "bitmap.h"
#include "util.h"
#include <algorithm>
#include <filesystem>
#include <mutex>
#include <thread>

SDL_Window* window = nullptr;
SDL_Surface* screen = nullptr;
std::mutex sdlLock, eventLock;
std::thread sdlUIThread;
std::vector<SDL_Event> savedEvents;
volatile static bool exitRequested = false;

/// try to create a frame window with the given dimensions
bool initGraphics(int frameWidth, int frameHeight)
{
	if (SDL_Init(SDL_INIT_VIDEO) < 0) {
		printf("Cannot initialize SDL: %s\n", SDL_GetError());
		return false;
	}
	window = SDL_CreateWindow("heX-Ray", SDL_WINDOWPOS_UNDEFINED,
							  SDL_WINDOWPOS_UNDEFINED, frameWidth,
						      frameHeight, SDL_WINDOW_SHOWN);
	if (!window) {
		printf("Window could not be created: %s\n", SDL_GetError());
		return false;
	}
	//SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_SOFTWARE);
	screen = SDL_GetWindowSurface(window);
	if (!screen) {
		printf("Cannot set video mode %dx%d - %s\n", frameWidth, frameHeight, SDL_GetError());
		return false;
	}
	void handleUIThread(void);
	sdlUIThread = std::thread(handleUIThread);
	return true;
}

/// closes SDL graphics
void closeGraphics(void)
{
	if (sdlUIThread.joinable()) sdlUIThread.join();
	SDL_Quit();
}

/// displays a VFB (virtual frame buffer) to the real framebuffer, with the necessary color clipping
void displayVFB(Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE])
{
	int rs = screen->format->Rshift;
	int gs = screen->format->Gshift;
	int bs = screen->format->Bshift;
	for (int y = 0; y < screen->h; y++) {
		Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
		for (int x = 0; x < screen->w; x++)
			row[x] = vfb[y][x].toRGB32(rs, gs, bs);
	}
	showUpdatedFullscreen();
}

// find an unused file name like 'hexray_0005.bmp'
static void findUnusedFN(char fn[], const char* suffix)
{
	int idx = 0;
	while (idx < 10000) {
		sprintf(fn, "hexray_%04d.%s", idx, suffix);
		if (!std::filesystem::exists(fn)) return;
		idx++;
	}
}

bool takeScreenshot(const char* filename)
{
	extern Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE]; // from main.cpp

	Bitmap bmp;
	bmp.generateEmptyImage(frameWidth(), frameHeight());
	for (int y = 0; y < frameHeight(); y++)
		for (int x = 0; x < frameWidth(); x++)
			bmp.setPixel(x, y, vfb[y][x]);
	bool res = bmp.saveImage(filename);
	if (res) printf("Saved a screenshot as `%s'\n", filename);
	else printf("Failed to take a screenshot\n");
	return res;
}

bool takeScreenshotAuto(Bitmap::OutputFormat fmt)
{
	char fn[256];
	findUnusedFN(fn, fmt == Bitmap::outputFormat_BMP ? "bmp" : "exr");
	return takeScreenshot(fn);
}

static bool isWindowRedrawEvent(const SDL_Event& ev)
{
	const Uint32 WINDOW_DAMAGED_EVENTS[] = {
		SDL_WINDOWEVENT_SHOWN, SDL_WINDOWEVENT_EXPOSED, SDL_WINDOWEVENT_MAXIMIZED, SDL_WINDOWEVENT_RESTORED
	};
	for (auto& id: WINDOW_DAMAGED_EVENTS) if (ev.window.event == id) return true;
	return false;
}

static bool handleSystemEvent(const SDL_Event& ev)
{
	switch (ev.type) {
		case SDL_QUIT:
			exitRequested = true;
			return true;
		case SDL_WINDOWEVENT:
		{
			if (isWindowRedrawEvent(ev)) showUpdatedFullscreen();
			return true;
		}
		case SDL_KEYDOWN:
		{
			switch (ev.key.keysym.sym) {
				case SDLK_ESCAPE:
					exitRequested = true;
					return true;
				case SDLK_F12:
					takeScreenshotAuto(Bitmap::outputFormat_BMP);
					return true;
			}
			break;
		}
	}
	return false;
}

void getSDLInputs(const Uint8*& keystate, int& mouseDeltaX, int& mouseDeltaY, std::vector<SDL_Event>& events)
{
	events.clear();
	eventLock.lock();
	events.swap(savedEvents);
	keystate = SDL_GetKeyboardState(nullptr);
	SDL_GetRelativeMouseState(&mouseDeltaX, &mouseDeltaY);
	eventLock.unlock();
}

// this thread is automatically run by initGraphics and supports early exit of the program
// by pressing the ESC key or closing the window using the "X" button. It also handles window
// redraw events, the F12 key (for screenshots), and saves all other events ("non-system") to
// a queue (savedEvents). The queue can be retrieved later (getSavedEvents)
void handleUIThread()
{
	SDL_Event ev;
	while (!exitRequested && SDL_WaitEvent(&ev)) {
		if (!handleSystemEvent(ev)) {
			eventLock.lock();
			savedEvents.push_back(ev);
			eventLock.unlock();
		}
	}
}

/// waits the user to indicate he/she wants to close the application (by either clicking on the "X" of the window,
/// or by pressing ESC). Since the actual event handling is done by the UI thread, we just wait for it to exit
void waitForUserExit(void)
{
	sdlUIThread.join(); // wait for the message loop to close
}

/// checks if the user indicated he/she wants to close the application (by either clicking on the "X" of the window,
/// or by pressing ESC)
bool checkForUserExit(void)
{
	return exitRequested;
}

/// returns the frame width
int frameWidth(void)
{
	if (screen) return screen->w;
	return 0;
}

/// returns the frame height
int frameHeight(void)
{
	if (screen) return screen->h;
	return 0;
}

void Rect::clip(int W, int H)
{
	x1 = std::min(x1, W);
	y1 = std::min(y1, H);
	w = std::max(0, x1 - x0);
	h = std::max(0, y1 - y0);
}

std::vector<Rect> getBucketsList(int bucketSize)
{
	std::vector<Rect> res;
	int BW = (frameWidth() - 1) / bucketSize + 1;
	int BH = (frameHeight() - 1) / bucketSize + 1;
	for (int y = 0; y < BH; y++) {
		for (int x = 0; x < BW; x++)
			res.push_back(Rect(x * bucketSize, y * bucketSize, (x + 1) * bucketSize, (y + 1) * bucketSize));
		if (y % 2) std::reverse(res.end() - BW, res.end()); // make the odd rows run right-to-left, not left-to-right
	}
	// clip the edge buckets to the frame dimensions:
	for (int i = 0; i < (int) res.size(); i++)
		res[i].clip(frameWidth(), frameHeight());
	return res;
}

bool drawRect(Rect r, const Color& c)
{
	r.clip(frameWidth(), frameHeight());

	int rs = screen->format->Rshift;
	int gs = screen->format->Gshift;
	int bs = screen->format->Bshift;

	Uint32 clr = c.toRGB32(rs, gs, bs);
	for (int y = r.y0; y < r.y1; y++) {
		Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
		for (int x = r.x0; x < r.x1; x++)
			row[x] = clr;
	}
	return true;
}

void showUpdatedFullscreen()
{
	sdlLock.lock();
	SDL_UpdateWindowSurface(window);
	sdlLock.unlock();
}

void showUpdated(Rect r)
{
	SDL_Rect sdlr = { r.x0, r.y0, r.w, r.h };
	sdlLock.lock();
	SDL_UpdateWindowSurfaceRects(window, &sdlr, 1);
	sdlLock.unlock();
}

bool displayVFBRect(Rect r, Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE])
{
	r.clip(frameWidth(), frameHeight());
	int rs = screen->format->Rshift;
	int gs = screen->format->Gshift;
	int bs = screen->format->Bshift;
	for (int y = r.y0; y < r.y1; y++) {
		Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
		for (int x = r.x0; x < r.x1; x++)
			row[x] = vfb[y][x].toRGB32(rs, gs, bs);
	}
	showUpdated(r);
	return true;
}

bool markRegion(Rect r, const Color& bracketColor)
{
	r.clip(frameWidth(), frameHeight());
	const int L = 8;
	if (r.w < L+3 || r.h < L+3) return true; // region is too small to be marked
	const Uint32 BRACKET_COLOR = bracketColor.toRGB32();
	const Uint32 OUTLINE_COLOR = Color(1, 1, 1).toRGB32();
	#define DRAW_ONE(x, y, color) \
		((Uint32*) (((Uint8*) screen->pixels) + ((r.y0 + (y)) * screen->pitch)))[r.x0 + (x)] = color
	#define DRAW(x, y, color) \
		DRAW_ONE(x, y, color); \
		DRAW_ONE(y, x, color); \
		DRAW_ONE(r.w - 1 - (x), y, color); \
		DRAW_ONE(r.w - 1 - (y), x, color); \
		DRAW_ONE(x, r.h - 1 - (y), color); \
		DRAW_ONE(y, r.h - 1 - (x), color); \
		DRAW_ONE(r.w - 1 - (x), r.h - 1 - (y), color); \
		DRAW_ONE(r.w - 1 - (y), r.h - 1 - (x), color)

	for (int i = 1; i <= L; i++) {
		DRAW(i, 0, OUTLINE_COLOR);
	}
	DRAW(1, 1, OUTLINE_COLOR);
	DRAW(L + 1, 1, OUTLINE_COLOR);
	for (int i = 0; i <= L; i++) {
		DRAW(i, 2, OUTLINE_COLOR);
	}
	for  (int i = 2; i <= L; i++) {
		DRAW(i, 1, BRACKET_COLOR);
	}

	showUpdated(r);

	return true;
}

/// displays pixels, set to true in the given array in yellow on the screen
void markAApixels(bool needsAA[VFB_MAX_SIZE][VFB_MAX_SIZE])
{
	Uint32 YELLOW = Color(1, 1, 0).toRGB32(screen->format->Rshift, screen->format->Gshift, screen->format->Bshift);
	for (int y = 0; y < screen->h; y++) {
		Uint32 *row = (Uint32*) ((Uint8*) screen->pixels + y * screen->pitch);
		for (int x = 0; x < screen->w; x++)
			if (needsAA[y][x]) row[x] = YELLOW;
	}
	showUpdatedFullscreen();
}

unsigned char SRGB_COMPRESS_CACHE[4097];

void Color::init_sRGB_cache(void)
{
	// precache the results of convertTo8bit_sRGB, in order to avoid the costly pow()
	// in it and use a lookup table instead, see Color::convertTo8bit_sRGB_cached().
	for (int i = 0; i <= 4096; i++)
		SRGB_COMPRESS_CACHE[i] = (unsigned char) convertTo8bit_sRGB(i / 4096.0f);
}

unsigned convertTo8bit_sRGB_cached(float x)
{
	if (x <= 0) return 0;
	if (x >= 1) return 255;
	return SRGB_COMPRESS_CACHE[int(x * 4096.0f)];
}
