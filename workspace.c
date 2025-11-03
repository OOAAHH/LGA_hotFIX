#include <stdlib.h>
#include "tools.h"

int lga_ws_init(lga_workspace **pws)
{
  if (!pws) return -1;
  if (*pws) return 0; /* already allocated */

  lga_workspace *ws = (lga_workspace*)calloc(1, sizeof(*ws));
  if (!ws) return -1;

  size_t n = (size_t)MAXRES + 1;
  size_t n2 = 2 * n;
  size_t n3 = 3 * n;

  /* Allocate required buffers */
  ws->calphas2 = (long*)calloc(n2, sizeof(long));
  ws->l_atoms2 = (long*)calloc(n2, sizeof(long));
  ws->out2     = (long*)calloc(n2, sizeof(long));
  ws->ind3     = (long*)calloc(3*n, sizeof(long));
  ws->tmp      = (long*)calloc(n, sizeof(long));
  ws->rwa1     = (long*)calloc(n, sizeof(long));
  ws->rwa2     = (long*)calloc(n, sizeof(long));
  ws->xa3      = (float*)calloc(n3, sizeof(float));
  ws->xb3      = (float*)calloc(n3, sizeof(float));
  ws->sup      = (long*)calloc(n2, sizeof(long));
  ws->sup1     = (long*)calloc(n2, sizeof(long));
  ws->sup2     = (long*)calloc(n2, sizeof(long));
  ws->mol      = (check_mol2*)calloc(2, sizeof(check_mol2));
  ws->straline = (stext*)calloc(200, sizeof(stext));

  if (!ws->calphas2 || !ws->l_atoms2 || !ws->out2 || !ws->ind3 || !ws->tmp ||
      !ws->rwa1 || !ws->rwa2 || !ws->xa3 || !ws->xb3 || !ws->sup ||
      !ws->sup1 || !ws->sup2 || !ws->mol || !ws->straline) {
    lga_ws_free(&ws);
    return -1;
  }

  *pws = ws;
  return 0;
}

void lga_ws_free(lga_workspace **pws)
{
  if (!pws || !*pws) return;
  lga_workspace *ws = *pws;
  free(ws->calphas2);
  free(ws->l_atoms2);
  free(ws->out2);
  free(ws->ind3);
  free(ws->tmp);
  free(ws->rwa1);
  free(ws->rwa2);
  free(ws->xa3);
  free(ws->xb3);
  free(ws->sup);
  free(ws->sup1);
  free(ws->sup2);
  free(ws->mol);
  free(ws->straline);
  free(ws);
  *pws = NULL;
}
