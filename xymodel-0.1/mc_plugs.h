struct plugin_stats
{
	size_t len;
	void (*get_data)(double *data, double (*lattice)[SIZE], double T);
	void (*get_x)(double *x);
	double *sumData;
	double *sumData2;
};

static struct plugin_stats pstat;

void update_statistic(const struct plugin_stats *stat, double (*lattice)[SIZE], double T)
{
	double data[stat->len];

	stat->get_data(data, lattice, T);

	size_t i;
	for (i=0; i<stat->len; i++)
	{
		stat->sumData[i] += data[i];
		stat->sumData2[i] += pow(data[i],2);
	}
}

void update_stat(double (*lattice)[SIZE], double T)
{
	update_statistic(&pstat, lattice, T);
}

void output_stat(FILE *out, struct plugin_stats *stat, size_t iterations)
{
	double x[stat->len];
	if (stat->get_x != NULL)
		stat->get_x(x);

	size_t i;
	for (i=0; i<stat->len; i++)
	{
		double sigma2 = 1.0/(iterations-1.0)
			* (stat->sumData2[i] - pow(stat->sumData[i],2)/iterations);

		if (stat->len > 4)
		{
			if (stat->get_x != NULL)
				fprintf(out, "%G\t", x[i]);
			else
				fprintf(out, "%d\t", i);
		}

		fprintf(out, "%G\t%G\n", stat->sumData[i]/iterations,
				sqrt(sigma2));
	}
}

void init_stat(char* lib)
{
	void *handle = dlopen(lib, RTLD_NOW);
	if (NULL == handle)
		error(1, 0, "error plug-in: %s", dlerror());

	dlerror(); /* clear errors */
	size_t (*get_length)() = dlsym(handle, "get_length");
	const char *err_msg;
	if ((err_msg = dlerror()) != NULL)
		error(1, 0, "get_length not found in %s: %s\n", lib, err_msg);

	pstat.len = get_length();
	pstat.sumData = (double *) malloc(pstat.len * sizeof(double));
	pstat.sumData2 = (double *) malloc(pstat.len * sizeof(double));

	size_t i;
	for (i=0; i<pstat.len; i++)
	{
		pstat.sumData[i] = 0;
		pstat.sumData2[i] = 0;
	}


	dlerror(); /* clear errors */
	pstat.get_data = dlsym(handle, "get_data");
	if ((err_msg = dlerror()) != NULL)
		error(1, 0, "get_data not found in %s: %s\n", lib, err_msg);

	dlerror(); /* clear errors */
	pstat.get_x = dlsym(handle, "get_x");
	if (dlerror() != NULL)
		pstat.get_x = NULL;
}
