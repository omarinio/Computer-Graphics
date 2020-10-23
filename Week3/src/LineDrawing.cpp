#include <CanvasTriangle.h>
#include <CanvasPoint.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <TextureMap.h>
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

	result.push_back(temp);

	for (int i = 0; i < numVals-1; i++) {
		temp = temp + step;
		result.push_back(temp);
	}
	return result;
}

std::vector<CanvasPoint> interpolatePoints(CanvasPoint start, CanvasPoint end, int steps){
    std::vector<CanvasPoint> points;
    for (int i=0; i<steps; i++){
        CanvasPoint p;
        if (steps == 1) steps = 2;
        p.x = start.x+((end.x-start.x)*i/(steps-1));
        p.y = start.y+((end.y-start.y)*i/(steps-1));
        points.push_back(p);
    }
    return points;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numVals) {
	std::vector<glm::vec3> result;
	float xStep = (to.x - from.x)/(numVals - 1);
	float yStep = (to.y - from.y)/(numVals - 1);
	float zStep = (to.z - from.z)/(numVals - 1);

	float xTemp = from.x;
	float yTemp = from.y;
	float zTemp = from.z;

	glm::vec3 first(xTemp, yTemp, zTemp);
	result.push_back(first);

	for (int i = 0; i < numVals-1; i++) {
		xTemp = xTemp + xStep;
		yTemp = yTemp + yStep;
		zTemp = zTemp + zStep;
		glm::vec3 temp(xTemp, yTemp, zTemp); 
		result.push_back(temp);
	}

	return result;
}

void drawLine(DrawingWindow &window, CanvasPoint from, CanvasPoint to, Colour colour) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numSteps;
	float yStepSize = yDiff/numSteps;
	for (float i = 0.0; i < numSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		uint32_t set = (255 << 24) + (colour.red << 16) + (colour.green << 8) + colour.blue;
		window.setPixelColour(round(x), round(y), set);
	}
}

void drawTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	drawLine(window, triangle[0], triangle[1], colour);
	drawLine(window, triangle[1], triangle[2], colour);
	drawLine(window, triangle[2], triangle[0], colour);
}

void fillHalf(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	CanvasPoint top = triangle[0];
	CanvasPoint mid = triangle[1];
	CanvasPoint bot = triangle[2];

	float x1Diff = mid.x - top.x;
	float y1Diff = mid.y - top.y;

	float x2Diff = bot.x - top.x;
	float y2Diff = bot.y - top.y;

	float steps1 = std::max(abs(x1Diff), abs(y1Diff));
	float steps2 = std::max(abs(x2Diff), abs(y2Diff));

	float x1Step = x1Diff/steps1;
	float x2Step = x2Diff/steps2;

	float y1Step = y1Diff/steps1;
	float y2Step = y2Diff/steps2;

	for (float i = 0.0; i < steps1; i++) {
		for (float j = 0.0; j < steps2; j++) {
			float x1 = top.x + (x1Step * i);
			float y1 = top.y + (y1Step * i);

			float x2 = top.x + (x2Step * j);
			float y2 = top.y + (y2Step * j);

			CanvasPoint point1(round(x1), round(y1));
			CanvasPoint point2(round(x2), round(y2));
			
			drawLine(window, point1, point2, colour);
		}
	}
}

void drawFilledTriangle(DrawingWindow &window, CanvasTriangle triangle, Colour colour) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	// getting flat coords

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	CanvasTriangle triangle1(top, mid, split);
	CanvasTriangle triangle2(mid, split, bot);

	fillHalf(window, triangle1, colour);
	fillHalf(window, triangle2, colour);
	drawTriangle(window, triangle, Colour(255,255,255));

}

void textureFill(DrawingWindow &window, CanvasTriangle triangle, TextureMap texture) {
	CanvasPoint top = triangle.vertices[0];
	CanvasPoint mid = triangle.vertices[1];
	CanvasPoint bot = triangle.vertices[2];

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	if (mid.y < top.y) {
		std::swap(mid.y, top.y);
		std::swap(mid.x, top.x);
	}

	if (bot.y < mid.y) {
		std::swap(bot.y, mid.y);
		std::swap(bot.x, mid.x);
	}

	CanvasPoint split;
	split.y = mid.y;

	split.x = round(top.x + ((mid.y - top.y)/(bot.y-top.y)) * (bot.x-top.x));

	split.texturePoint = TexturePoint(round(top.texturePoint.x + ((mid.texturePoint.y - top.texturePoint.y)/(bot.texturePoint.y-top.texturePoint.y)) * (bot.texturePoint.x-top.texturePoint.x)), 
		round(top.texturePoint.y + ((mid.texturePoint.y - top.texturePoint.y)/(bot.texturePoint.y-top.texturePoint.y)) * (bot.texturePoint.y-top.texturePoint.y)));


	// TOP TRIANGLE ---------------------------------------------------------------------------------------------------------------------------------------------------
	int height = mid.y - top.y;

	std::vector<CanvasPoint> starts = interpolatePoints(top, mid, height);
	std::vector<CanvasPoint> ends = interpolatePoints(top, split, height);
	std::vector<CanvasPoint> textStarts = interpolatePoints(CanvasPoint(top.x, top.y), CanvasPoint(mid.x, mid.y), height);
	std::vector<CanvasPoint> textEnds = interpolatePoints(CanvasPoint(top.x, top.y), CanvasPoint(split.x, split.y), height);

	for (int y = 0; y < height; y++) {
		int width = ends.at(y).x - starts.at(y).x;

		std::vector<CanvasPoint> line = interpolatePoints(textStarts.at(y), textEnds.at(y), width);
		for (int x=0; x < width; x++) {
			int x_coord = int(line.at(x).x);
			int y_coord = int(line.at(x).y);
			uint32_t col = texture.pixels.at((y_coord*texture.width) + x_coord);
			window.setPixelColour(starts.at(y).x+x, starts.at(y).y, col);
		}
	}

	// BOTTOM TRIANGLE ------------------------------------------------------------------------------------------------------------------------------------------------
	int height2 = bot.y - mid.y;

	std::vector<CanvasPoint> starts2 = interpolatePoints(mid, bot, height2);
	std::vector<CanvasPoint> ends2 = interpolatePoints(split, bot, height2);
	std::vector<CanvasPoint> textStarts2 = interpolatePoints(CanvasPoint(mid.x, mid.y), CanvasPoint(bot.x, bot.y), height2);
	std::vector<CanvasPoint> textEnds2 = interpolatePoints(CanvasPoint(split.x, split.y), CanvasPoint(bot.x, bot.y), height2);

	for (int y = 0; y < height2; y++) {
		int width = ends2.at(y).x - starts2.at(y).x;

		std::vector<CanvasPoint> line = interpolatePoints(textStarts2.at(y), textEnds2.at(y), width);
		for (int x=0; x < width; x++) {
			int x_coord = int(line.at(x).x);
			int y_coord = int(line.at(x).y);
			uint32_t col = texture.pixels.at((y_coord*texture.width) + x_coord);
			window.setPixelColour(starts2.at(y).x+x, starts2.at(y).y-1, col);
		}
	}

}

void draw(DrawingWindow &window) {
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = rand() % 256;
			float green = 0.0;
			float blue = 0.0;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
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
		else if (event.key.keysym.sym == SDLK_u) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand()%window.width, rand()%window.height), 
				CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));

			Colour colour = Colour(rand()%256, rand()%256, rand()%256);
			drawTriangle(window, triangle, colour); 
		}
		else if (event.key.keysym.sym == SDLK_f) {
			CanvasTriangle triangle = CanvasTriangle(CanvasPoint(rand()%window.width, rand()%window.height), 
				CanvasPoint(rand()%window.width, rand()%window.height), CanvasPoint(rand()%window.width, rand()%window.height));

			Colour colour = Colour(rand()%256, rand()%256, rand()%256);
			drawFilledTriangle(window, triangle, colour); 
		}
		else if (event.key.keysym.sym == SDLK_t) {
			TextureMap texture = TextureMap("texture.ppm");
			CanvasPoint p1(160, 10);
			CanvasPoint p2(300,230);
			CanvasPoint p3(10, 150);

			p1.texturePoint = TexturePoint(195, 5);
			p2.texturePoint = TexturePoint(395, 380);
			p3.texturePoint = TexturePoint(65, 330); 

			CanvasTriangle triangle(p1, p2, p3);

			textureFill(window, triangle, texture);

		}
		
	} else if (event.type == SDL_MOUSEBUTTONDOWN) window.savePPM("output.ppm");
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;

	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		update(window);
		draw(window);

		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}

