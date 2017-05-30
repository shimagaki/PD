#include <unistd.h>
#include <stdlib.h>
#include <argp.h>
#include <math.h>
#include <string.h>

#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h> 
#endif

#include "lattice.h"
#include "vortex.h"

const char *argp_program_version = "glvec 0.1";
static char doc[] = 
"Plot spins. Use [+] and [-] keys to zoom; move with hjkl. "
"Press v to toggle highlighting of vortices and e to toggle between vectors and energy.\n"
"Quit with [Esc]"
;

static char args_doc[] = "FILE";

struct arguments
{
	char* shm_file;
	int update_interval; //msecs
};

static struct argp_option options[] = {
	{"update-interval", 'u', "SEC", 0,
		"Update every SEC seconds. (default 1.0)" },
	{"draw-vortices",   'v',  0,    0,
		"Show vortices. blue: clockwise, red: counterclockwise"},
	{"energy",          'e',  0,    0, "Draw energy of spin"},
	{"block",            'b',  0,    0, "Draw colored blocks"},
	{ 0 }
};

static struct d_options
{
	enum  {vectors, energy, block} class;
	int vortices;
} draw_options = {vectors, 0};

static error_t parse_opt(int key, char *arg, struct argp_state *state)
{
	struct arguments *arguments = state->input;

	switch (key)
	{
		case 'u':
			errno = 0;
			arguments->update_interval = strtod(arg, NULL) * 1E3;
			if (errno != 0)
			{
				perror("update_interval");
				argp_usage(state);
			}
			break;
		case 'v':
			draw_options.vortices = 1;
			break;
		case 'e':
			draw_options.class = energy;
			break;
		case 'b':
			draw_options.class = block;
			break;
		case ARGP_KEY_ARG:
			if (state->arg_num > 0)
				argp_usage(state);
			arguments->shm_file = arg;
			break;
		case ARGP_KEY_END:
			if (state->arg_num < 1)
				argp_usage(state);
			break;
	}

	return 0;
}

static struct argp argp = {options, parse_opt, args_doc, doc};

int window; 

int mid_x;
int mid_y;
int zoom = 1;

int is_visible = GLUT_VISIBLE;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

#define X_LEFT (mid_x - SIZE/2/zoom)
#define X_RIGHT (mid_x + SIZE/2/zoom)

#define Y_LOWER (mid_y - SIZE/2/zoom)
#define Y_UPPER (mid_y + SIZE/2/zoom)

GLuint CIRCLE_LIST;
GLuint ARROW_LIST;

void *lattice_shm;

void init_display_lists()
{
	CIRCLE_LIST = glGenLists(2);
	ARROW_LIST = CIRCLE_LIST + 1;

	const size_t shapes = 180;

	glNewList(CIRCLE_LIST, GL_COMPILE);
	glBegin(GL_POLYGON);

	GLint i;
	for(i=0;i<shapes;i++)
		glVertex2f(cos(i*2*M_PI/shapes), sin(i*2*M_PI/shapes));

	glEnd();
	glEndList();


	glNewList(ARROW_LIST, GL_COMPILE);
	glBegin(GL_TRIANGLES);
	glVertex2d( 0.0f, 0.4f);
	glVertex2d( 0.1f, 0.2f);
	glVertex2d(-0.1f, 0.2f);
	glEnd();

	glBegin(GL_LINES);
	glVertex2d(0, 0.2);
	glVertex2d(0, -0.4);
	glEnd();
	glEndList();
}

void init_GL()
{
	init_display_lists();

	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glClearDepth(1.0);
	glShadeModel(GL_SMOOTH);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);                                                    
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	mid_x = SIZE/2;
	mid_y = SIZE/2;

	glOrtho(0, SIZE, 0, SIZE,  1, -1);
	glMatrixMode(GL_MODELVIEW);
}

//TODO please no copy&paste
double getH(double (*lattice)[SIZE], int x, int y)
{
	const double spin = lattice[x][y];
	double H;
	H = cos(spin - lattice[(x-1)&PBC][y]);
	H += cos(spin - lattice[(x+1)&PBC][y]);
	H += cos(spin - lattice[x][(y+1)&PBC]);
	H += cos(spin - lattice[x][(y-1)&PBC]);

	return H/4.0;
}

float interpolate(double alpha, double rotate)
{
	float x =  (alpha + rotate) / (2 * M_PI);

	if (x >= 1) // x [0,1)
		x -= (floor(x));

	x = 2 - 3 *fabsf(2*x - 1);

	if (x > 1)
		return 1;
	else if (x < 0)
		return 0;

	return x;
}

void draw_block(double (*lattice)[SIZE], int draw_energy)
{
	glDisable(GL_POLYGON_SMOOTH);

	size_t x;
	for (x=X_LEFT; x < X_RIGHT; x++) 
	{
		size_t y;
		for (y=Y_LOWER; y< Y_UPPER; y++) 
		{
			glBegin(GL_QUADS);
			if (draw_energy)
			{
				double color = getH(lattice, x, y);
				glColor3f(1.0f-color, 0.0f, color);
			}
			else
			{
				float r = interpolate(lattice[x][y], M_PI);
				float g = interpolate(lattice[x][y], 5.f/3*M_PI);
				float b = interpolate(lattice[x][y], 1.f/3*M_PI);
				glColor3f(r,g,b);
			}
			glVertex2d(x, y);
			glVertex2d(x, y + 1.0f);
			glVertex2d(x + 1.0f, y + 1.0f);
			glVertex2d(x + 1.0f, y);
			glEnd();
		}
	}
	glEnable(GL_POLYGON_SMOOTH);
}

void draw_vectors(double (*lattice)[SIZE])
{
	size_t x;
	for (x=X_LEFT; x < X_RIGHT; x++) 
	{
		size_t y;
		for (y=Y_LOWER; y< Y_UPPER; y++) 
		{
			double angle = lattice[x][y]/M_PI*180;

			glTranslatef(x+0.5f, y+0.5f, 0); 
			glRotatef(angle, 0.0f, 0.0f, 1.0f);

			float color_angle = fabsf(1.0f - lattice[x][y]/(M_PI));
			glColor3f(1.0f - color_angle, 0.0f , color_angle);  

			glCallList(ARROW_LIST);

			glLoadIdentity();
		}
	}

}

void draw_vortices(double (*lattice)[SIZE])
{
	struct point p = next_vortex(lattice);
	while( p.type != -1)
	{
		if (CLOCKWISE == p.type)
			glColor4f(0, 0, 1, 0.4);
		else
			glColor4f(1, 0, 0, 0.4);

		p.x += 1.5;
		p.y += 1.5;

		glTranslatef(p.x,p.y,0);
		glCallList(CIRCLE_LIST);
		glLoadIdentity();

		p = next_vortex(lattice);
	}

}

void draw()
{
	glClear(GL_COLOR_BUFFER_BIT);	
	glLoadIdentity();

	double (*lattice)[SIZE] = lattice_shm;

	switch(draw_options.class)
	{

		case energy:
			draw_block(lattice, 1);
			break;

		case block:
			draw_block(lattice, 0);
			break;

		case vectors:
			draw_vectors(lattice);
			break;
	}

	if (draw_options.vortices)
		draw_vortices(lattice);

	glutSwapBuffers();
}

static void update_zoom()
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glOrtho(X_LEFT, X_RIGHT, Y_LOWER, Y_UPPER, 1, -1);
	glMatrixMode(GL_MODELVIEW);

	if (GLUT_VISIBLE == is_visible)
		glutPostRedisplay();
}

void key_pressed(unsigned char key, int x, int y) 
{
	/* avoid thrashing this procedure */
	usleep(100);

	switch (key)
	{
		case '\e':	/* Esc  (Quit) */
		case 'q':
			glutDestroyWindow(window); 
			exit(0);                   


		/* zoom */
		case '+':
			if (SIZE/2/zoom >= 2)
			{
				zoom *= 2;
				update_zoom();
			}
			break;
		case '-':
			if (zoom <= 1)
				zoom = 1;
			else
			{
				zoom /= 2;
				if (mid_x < SIZE/2/zoom)
					mid_x = SIZE/2/zoom;
				else if (mid_x > SIZE - SIZE/2/zoom )
					mid_x = SIZE - SIZE/2/zoom;

				if (mid_y < SIZE/2/zoom)
					mid_y = SIZE/2/zoom;
				else if (mid_y > SIZE - SIZE/2/zoom)
					mid_y = SIZE - SIZE/2/zoom;

				update_zoom();
			}
			break;
		case 'h':	/* Left */
			mid_x = MAX(mid_x - 1, SIZE/2/zoom);
			update_zoom();
			break;
		case 'l':	/* Right */
			mid_x = MIN(mid_x + 1, SIZE-SIZE/2/zoom);
			update_zoom();
			break;
		case 'j':	/* Down */
			mid_y = MAX(mid_y - 1, SIZE/2/zoom);
			update_zoom();
			break;
		case 'k':	/* Up */
			mid_y = MIN(mid_y + 1, SIZE-SIZE/2/zoom);
			update_zoom();
			break;


		/* draw options */
		case 'v':
			draw_options.vortices = !draw_options.vortices & 1;
			glutPostRedisplay();
			break;
		case 'e':
			draw_options.class = energy;
			glutPostRedisplay();
			break;
		case 'b':
			draw_options.class = block;
			glutPostRedisplay();
			break;
		case 'a':
			draw_options.class = vectors;
			glutPostRedisplay();
			break;
	}
}


void update_visible(int state)
{
	is_visible = state;
}

void redraw(int update_interval)
{
	if (GLUT_VISIBLE == is_visible)
		glutPostRedisplay();
	glutTimerFunc(update_interval, &redraw, update_interval);
}

int main(int argc, char *argv[]) 
{  
	glutInit(&argc, argv);  

	struct arguments arguments;
	arguments.shm_file = NULL;
	arguments.update_interval = 1E3;

	argp_parse(&argp, argc, argv, 0, 0, &arguments);

	void * shmaddr = get_shm_addr(arguments.shm_file);

	readXY(shmaddr);
	lattice_shm = get_shm_lattice(shmaddr);


	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);  

	glutInitWindowSize(500, 500);  
	
	const char title_format[] = "XY model (%s)";
	char title[sizeof(title_format) + strlen(arguments.shm_file) + 1];
	snprintf(title, sizeof(title), title_format, arguments.shm_file);
	window = glutCreateWindow(title);  

	glutVisibilityFunc(&update_visible);
	glutDisplayFunc(&draw);  
	glutKeyboardFunc(&key_pressed);

	init_GL();

	glutTimerFunc(arguments.update_interval, &redraw, arguments.update_interval);

	glutMainLoop();  

	exit(1); /* not reached */
}
