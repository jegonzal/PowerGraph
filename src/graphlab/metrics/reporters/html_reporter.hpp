
#ifndef GRAPHLAB_HTML_REPORTER
#define GRAPHLAB_HTML_REPORTER

#include <graphlab/metrics/imetrics_reporter.hpp>
 

/**
  * Simple metrics reporter that dumps metrics to HTML
  */

namespace graphlab {

    class html_reporter : public imetrics_reporter {
        private:
            html_reporter() {}
        
            std::string filename;
            FILE * f;
        public:
            
        
            html_reporter(std::string fname) : filename(fname) {
                 // Create new file
                 f = fopen(fname.c_str(), "w");
                 assert(f != NULL);
                 fprintf(f, "<html><head><title>GraphLab Metrics Report</title>");
                 fprintf(f, "<style>\n");
                 fprintf(f, "table {  border: 1px solid #999999; font: normal 80%%/140%% arial, helvetica, sans-serif; color: #555; background: #fff;}  td, th {border: 1px dotted #bbb; padding: .5em; width:100px}   ");
                 fprintf(f, "</style></head><body>");
            }
            
            ~html_reporter() {
                fprintf(f, "</body></html>");
                fclose(f);
            }
            
            virtual void do_report(std::string name, std::string ident, std::map<std::string, metrics_entry> & entries) {
                if (ident != name) {
                    fprintf(f, "<h3>%s:%s</h3>\n", name.c_str(), ident.c_str());
                } else {
                    fprintf(f, "<h3>%s</h3>\n", name.c_str());
                }
                
                
                // First write numeral, then timings, then string entries
                for(int round=0; round<3; round++) { 
                    std::map<std::string, metrics_entry>::iterator it;
                    int c = 0;
                    fprintf(f, "<!-- Round %d -->\n", round);
                    fprintf(f, "\n<p>");
                    for(it = entries.begin(); it != entries.end(); ++it) {
                        metrics_entry ent = it->second;
                        switch(ent.valtype) {
                            case INTEGER:
                                if (round == 0) {   
                                    if (c++ == 0)
                                        fprintf(f, "<table><tr><th>Key</th><th>Value</th><th>Count</th><th>Min</th><th>Max</th><th>Average</th></tr>");
                                        
                                    fprintf(f, "<tr><td>%s</td>\n",  it->first.c_str());
                                    
                                    fprintf(f, "<td>%ld</td>\n", (long int) ent.value);
                                    if (ent.count > 1) {
                                        fprintf(f, "<td>%ld</td>\n", (long int) ent.count);
                                        fprintf(f, "<td>%ld</td>\n", (long int) ent.minvalue);
                                        fprintf(f, "<td>%ld</td>\n", (long int) ent.maxvalue);
                                        fprintf(f, "<td>%.3lf</td>\n",  ent.cumvalue/ent.count);
                                     } else fprintf(f, "<td colspan=4>&nbsp;</td>");
                                     fprintf(f, "</tr>");
                                }
                                break;
                            case REAL:
                               if (round == 0) {   
                                    if (c++ == 0)
                                        fprintf(f, "<table><tr><th>Key</th><th>Value</th><th>Count</th><th>Min</th><th>Max</th><th>Average</th></tr>");
                               }
                            case TIME:
                                if (ent.valtype == TIME && round == 1) {
                                    if (c++ == 0) 
                                         fprintf(f, "<table><tr><th>Key</th><th>Value (sec)</th><th>Count</th><th>Min (sec)</th><th>Max (sec)</th><th>Average (sec)</th></tr>\n");
                                }
                                if ((round == 0 && ent.valtype == REAL)||(round == 1 && ent.valtype == TIME)) {
                                    fprintf(f, "<tr><td>%s</td>\n",  it->first.c_str());
                                    
                                    fprintf(f, "<td>%lf</td>\n",  ent.value);
                                    if (ent.count > 1) {
                                        fprintf(f, "<td>%ld</td>\n", (long int) ent.count);
                                        fprintf(f, "<td>%.3lf</td>\n",  ent.minvalue);
                                        fprintf(f, "<td>%.3lf</td>\n",  ent.maxvalue);
                                        fprintf(f, "<td>%.3lf</td>\n",  ent.cumvalue/ent.count);
                                    } else fprintf(f, "<td colspan=4>&nbsp;</td>");
                                    fprintf(f, "</tr>");
                                } 
                                break;
                            case STRING:
                                if (round == 2) {
                                    if (c++ == 0)
                                        fprintf(f, "<table><tr><th>Key</th><th>Value</th></tr>\n"); 
                                    fprintf(f, "<tr><td>%s</td><td width=400>%s</td>\n",  it->first.c_str(), ent.stringval.c_str());                                
                                    fprintf(f, "</tr>");
                                }
                                
                                break;
                        }
                    }
                    if (c>0) fprintf(f, "</table>");
                    fprintf(f, " </p>");
                }

                fflush(f);
            };
        
    };
    
};



#endif



