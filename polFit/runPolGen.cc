#include "general/ArgParser.h"

#include "polGen.C"

#include <string>

#if !(defined(__CINT__) or defined(__CLING__))
int main(int argc, char *argv[])
{
  ArgParser parser(argc, argv);
  const auto outfile = parser.getOptionVal<std::string>("--outfile");
  const auto lthsig = parser.getOptionVal<double>("--lthsig");
  const auto lthbkg = parser.getOptionVal<double>("--lthbkg");
  const auto lphsig = parser.getOptionVal<double>("--lphsig");
  const auto lphbkg = parser.getOptionVal<double>("--lphbkg");
  const auto ltpsig = parser.getOptionVal<double>("--ltpsig");
  const auto ltpbkg = parser.getOptionVal<double>("--ltpbkg");
  const auto nEvents = parser.getOptionVal<long>("--nEvents");

  polGen(outfile, lthsig, lthbkg, lphsig, lphbkg, ltpsig, ltpbkg, nEvents);

  return 0;
}

#endif
