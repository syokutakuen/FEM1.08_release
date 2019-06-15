/* 要素情報構造体の表示 */
void dump_eiinfo(FILE *fp, const ElemIntegralInfo *p)
{
	int i, j, k;
	fprintf (fp, "%-40s:%5d\n%-40s:%5d\n%-40s:%5d\n%-40s:%5d\n%-40s:%5d\n",
			"Data number", p->number,
			"Number of nodes", p->nnode,
			"Number of order", p->norder,
			"Number of total integral point count", p->ntint,
			"Number of dimension of shape function", p->ndim);
	fputs ("Shape function at integral points\n", fp);
	for (i = 0; i < p->ntint; i++) {
		double *dp;
		fprintf (fp, "\t%-40s:", "Coordinates at integral points");
		for (j = 0, dp = p->ipt + p->ndim * i; j < p->ndim; j++, dp++)
			fprintf (fp, "%15g", *dp);
		fputc ('\n', fp);
		fprintf (fp, "\t%-40s:", "Weight at integral points");
		fprintf (fp, "%15g", p->wt [i]);
		fputc ('\n', fp);
		if (p->N != NULL) {
			fputs ("\tShape function at integral points\n\t", fp);
			for (j = 0, dp = p->N + p->nnode * i; j < p->nnode; j++, dp++)
				fprintf (fp, "%15g", *dp);
			fputc ('\n', fp);
		}
		if (p->dN != NULL) {
			fputs ("\tDifference shape function at integral points\n\t", fp);
			for (j = 0, dp = p->dN + p->nnode * p->ndim * i; j < p->ndim; j++) {
			 for (k = 0; k < p->nnode; k++, dp++)
			 	fprintf (fp, "%15g", *dp);
			 fputc ('\n', fp);
			}
		}
	}
	fputc ('\n', fp);
}

void dump_einfo (FILE *fp, const ElementInfo *p)
{
	fprintf (fp, "%-50s:%d\n" "%-50s:%s\n" "%-50s:%d\n" "%-50s:%d\n"
		 "%-50s:%d\n" "%-50s:%d\n" "%-50s:%d\n\n",
		 "Data number", p->number,
		 "Data name", p->name,
		 "Number of vertex", p->nvert,
		 "Number of faces", p->nface,
		 "Reference to information number at body(main)", p->info1->number,
		 "Reference to information number at body(sub)", (p->info2 == NULL)? 0: p->info2->number,
		 "Reference to information number at face", p->info3->number);
}

void dump_all_einfo (FILE *fp, int flag)
{
	int i;
	if (!flag) return;
	fputs ("\nSummary of Element caluration information table\n", fp);
	for (i = 1; i < neiinfo; i++)
		dump_eiinfo (fp, eiinfo + i);
	fputs ("\nSummary of Element information table\n", fp);
	for (i = 1; i < neinfo; i++)
		dump_einfo (fp, einfo + i);
}

/* 要素の節点座標値をマトリックス表示する。 デバッグ用 */
void dump_X_mat (const Element *e, int flag)
{
	int j, k;
	int ndim = etinfo [prob].ndim, nnode = e->info->info1->nnode;
	if (!flag) return;
	for (j = 0; j < nnode; j++) {
		for (k = 0; k < ndim; k++)
			fprintf (stderr, "%15g", coord [e->conn [j] * mdim + k]);
		fputc ('\n', stderr);
	}
	fputc ('\n', stderr);
}

/* 2次元配列のダンプ */
void dump_mat2 (const double *p, int sz1, int sz2, int f)
{
	int i, j;
	if (!f) return;
	for (i = 0; i< sz1; i++) {
		for (j = 0; j < sz2; j++)
			fprintf (stderr, "%15g", *p++);
		fputc ('\n', stderr);
	}
	fputc ('\n', stderr);
}

/* 全体剛性のダンプ */
void dump_pline_info (int f)
{
	int i;
	K_line_info *p;

	if (!f) return;
	fprintf (stderr, "%5s%5s%6s%6s%6s%6s%6s\n",
		"node", "dof", "pos", "strat", "end", "size", "pos2");
	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++)
		fprintf (stderr, "%6d%6d%6d%6d" LLU(6) "\n",
				p->pos, p->start, p->end, p->size, p->p - sysk);
}

void dump_pline_data (int f)
{
	int i, j;
	K_line_info *p;

	if (!f) return;
	fputs ("\ntotal stiff matrix\n", stderr);
	for (i = 0, p = line_info; i < rank_line_info [KLDOF_NONE][1]; i++, p++) {
		for (j = 0; j < p->start; j++)
			fprintf (stderr, "%15c", ' ');
		for (j = p->start; j < p->end; j++)
			fprintf (stderr, "%15g", line_info [i].p [j]);
		fputc ('\n', stderr);
	}
	fputc ('\n', stderr);
}

void dump_pline (int f1, int f2)
{
	dump_pline_info (f1);
	dump_pline_data (f2);
}

void dump_rev_line_info (int f)
{
	int i, pos;
	if (!f) return;
	fputs ("dump_rev_line_info\n", stderr);
	for (i = 0; i < ntnode * mdof; i++) {
		pos = rev_line_info [i];
		fprintf (stderr, "%5d%5d%5d%5d", i, pos, pos / mdof, pos % mdof);
		pos = decode_rev_info [i];
		fprintf (stderr, "%5d%5d%5d\n", pos, pos / mdof, pos % mdof);
	}
	fputc ('\n', stderr);
}

void dump_bc (int f)
{
	int i;
	if (!f) return;
	fputs ("dump_bc\n", stderr);
	for (i = 0; i < rank_line_info [KLDOF_NONE][1]; i++)
		fprintf (stderr, "%5d%15lg\n", i, rhs_value [i]);
	fputc ('\n', stderr);
}

void dump_nodal_dof (int f)
{
	int i;
	if (!f) return;
	fputs ("dump_nodal_dof\n", stderr);
	for (i = 0; i < ntnode * mdof; i++) {
		if (i > 0 && i % 20 == 0) fputc ('\n', stderr);
		fprintf (stderr, "%5d", nodal_dof_list [i]);
	}
	fputc ('\n', stderr);
}

void dump_rank_info (int f)
{
	int kltype;
	if (!f) return;
	fputs ("dump_rank_info\n", stderr);
	for (kltype = KLDOF_NONE; kltype <= KLDOF_ORPHAN; kltype++)
		fprintf (stderr, "%5d%5d%5d\n", kltype, rank_line_info [kltype][0], rank_line_info [kltype][1]);
	fputc ('\n', stderr);
}


