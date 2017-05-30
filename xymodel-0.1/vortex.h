#ifndef FINDVORTEX_H
#define FINDVORTEX_H

struct point
{
	size_t x,y;
	enum {CLOCKWISE, COUNTERCLOCKWISE} type;
};

static struct point p = {0,0,-1};
void reset_next_vortex()
{
	p.x = 0;
	p.y = 0;
	p.type = -1;
}

/* Returns next vortex. point.type =-1 if none left,
 * next call returns first one. */
struct point next_vortex(double (*lattice)[SIZE])
{
	p.type = -1;
	size_t x = p.x;
	for (; x < SIZE; x++) 
	{
		size_t y;
		if (x == p.x)
			y = p.y+1;
		else
			y = 0;
		for (; y < SIZE; y++) {
			double last = lattice[x][y];	// init
			int countNeg = 0;
			int countPos = 0;
			size_t i;
			// go around clockwise
			for (i=1; i<5; i++) {
				int dx = (i/2)%2;
				int dy = (i/2+i%2)%2;
				double alpha = lattice[(x+dx)&PBC][(y+dy)&PBC];

				// get the angle between the vectors (+ is clockwise)
				double angle = alpha-last;
				if (angle < -M_PI) angle += 2*M_PI;
				else if (angle > M_PI) angle -= 2*M_PI;

				if (fabs(angle) < M_PI) {
					if (angle >= 0) countPos++;
					if (angle <= 0) countNeg++;
				}
				else break;

				last = alpha;	// prepare for next compare
			}
			p.x = x;
			p.y = y;

			if(4 == countPos && 0 == countNeg) {
				p.type = CLOCKWISE;
				return p;
			}
			else if (4 == countNeg && 0 == countPos) {
				p.type = COUNTERCLOCKWISE;
				return p;
			}
		}
	}
	 
	reset_next_vortex();
	
	return p;
}

#endif
