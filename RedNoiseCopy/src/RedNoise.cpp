#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>

#define WIDTH 10000
#define HEIGHT 10000

void draw(DrawingWindow &window) {
	window.clearPixels();
	for (size_t y = 0; y < window.height / 20; y++) {
		for (size_t x = 0; x < window.width / 20; x++) {
			float red = 0.0;
			float green = 0.0;
			float blue = 0.0;
			double brightness =  fmax(0.5, sin((rand() % 1000000000) / 1000000000.0));
			red = 255 * brightness;
			green = 255 * brightness;
			blue = 255 * brightness;
			
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			for (int i = 0; i < 20; i++) {
				for (int f = 0; f < 20; f++) {
					window.setPixelColour(x * 20 + i, y * 20 + f, colour);
				}
			}
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
