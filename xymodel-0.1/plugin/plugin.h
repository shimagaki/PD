#ifndef PLUGIN_H
#define PLUGIN_H

__attribute__ ((const))
inline size_t get_length();

void get_data(double *cor, void *l, double T);

__attribute__ ((unused))
void get_x(double *x);


__attribute__ ((visibility("hidden")))
int main(int argc, char *argv[]);

#endif
