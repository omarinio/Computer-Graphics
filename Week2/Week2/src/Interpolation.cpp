#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/vec3.hpp> // glm::vec3
#include <glm/gtx/string_cast.hpp>

#define WIDTH 320
#define HEIGHT 240

std::vector<float> interpolateSingleFloats(float from, float to, int numVals) {
	std::vector<float> result;
	float step = (to - from)/(numVals-1); 
	float temp = from; 
	for (int i = 0; i < numVals; i++) {
		result.push_back(temp);
		temp = temp + step;
	}
	return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numVals) {
	std::vector<glm::vec3> result;
	float xStep = (to.x - from.x)/(numVals - 1);
	float yStep = (to.y - from.y)/(numVals - 1);
	float zStep = (to.z - from.z)/(numVals - 1);

	float xTemp = from.x;
	float yTemp = from.y;
	float zTemp = from.z;

	for (int i = 0; i < numVals; i++) {
		glm::vec3 temp(xTemp, yTemp, zTemp); 
		result.push_back(temp);
		xTemp = xTemp + xStep;
		yTemp = yTemp + yStep;
		zTemp = zTemp + zStep;
	}

	return result;
}

void draw(DrawingWindow &window, std::vector<float> weight) {
	window.clearPixels();
	for (size_t x = 0; x < window.width; x++) {
		for (size_t y = 0; y < window.height; y++) {
			float red = weight.at(x);
			float green = weight.at(x);
			float blue = weight.at(x);
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}

void update(DrawingWindow &window) {
	// Function for performing animation (shifting artifacts or moving the camera)
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	std::vector<glm::vec3> result;
	glm::vec3 from(1, 4, 9.2);
	glm::vec3 to(4, 1, 9.8);

	result = interpolateThreeElementValues(from, to, 4);
	for(size_t i=0; i<result.size(); i++) std::cout << glm::to_string(result[i]) << " ";
	std::cout << std::endl;

	// DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	// SDL_Event event;
	// while (true) {
	// 	// We MUST poll for events - otherwise the window will freeze !
	// 	if (window.pollForInputEvents(event)) handleEvent(event, window);
	// 	update(window);
	// 	draw(window, result);
	// 	// Need to render the frame at the end, or nothing actually gets shown on the screen !
	// 	window.renderFrame();
	// }
}
