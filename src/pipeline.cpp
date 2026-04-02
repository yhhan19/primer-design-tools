#include "pipeline.hpp"

void PDRStage::run(PipelineContext& ctx) {
    std::cout << "Running PDR optimizer\n";
    srand(ctx.args.seed);

    RiskOptimizer ro(ctx.labels,
                     ctx.sequences, 
                     ctx.args.len_PDR, 
                     ctx.args.len_amp, 
                     ctx.args.len_amp);
    auto PDR = ro.search(0, ctx.args.u_max, ctx.args.u_min);
    ro.validate_PDR(PDR);
    risk_t score = ro.score(PDR, ALPHA);
    std::cout << "loss: " << score << std::endl;

    std::ofstream fout(ctx.args.output_file + "." + name());
    fout << display_PDR(PDR) << std::endl;
    fout.close();

    ctx.pdr_regions = PDR;
}

void PrimerSelStage::run(PipelineContext& ctx) {

}

void DimerStage::run(PipelineContext& ctx) {
    std::cout << "Running dimer optimization\n";
    srand(ctx.args.seed);
    KPartiteGraph g(ctx.args.input_file);
    auto solution = g.solve_fast(ctx.args.iter);
    std::cout << "loss: " << g.cost(solution) << std::endl;

    std::ofstream fout(ctx.args.output_file);
    for (index_t i : solution)
        fout << i << " ";
    fout << "\n";
    fout.close();

    ctx.dimer_solution = solution;
}

void OffTargetStage::run(PipelineContext& ctx) {
    
}
