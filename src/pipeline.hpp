#pragma once
#include "utility.hpp"
#include "automaton.hpp"
#include "risk_optimizer.hpp"
#include "graph.hpp"

class Stage {
public:
    virtual ~Stage() = default;
    virtual void run(PipelineContext& ctx) = 0;
    virtual std::string name() const = 0;
    virtual std::string short_name() const = 0;
};

class PDRStage       : public Stage {
public:
    std::string name() const { return "PDR"; }
    std::string short_name() const { return "pdr"; }
    void run(PipelineContext& ctx);
    ~PDRStage() {}
};

class PrimerSelStage : public Stage {
public:
    std::string name() const { return "PrimerSelection"; }
    std::string short_name() const { return "p3"; }
    void run(PipelineContext& ctx);
    ~PrimerSelStage() {}
};

class DimerStage     : public Stage {
public:
    std::string name() const { return "DimerMinimization"; }
    std::string short_name() const { return "dim"; }
    void run(PipelineContext& ctx);
    ~DimerStage() {}
};

class OffTargetStage : public Stage {
public:
    std::string name() const { return "OffTargetSearch"; }
    std::string short_name() const { return "offt"; }
    void run(PipelineContext& ctx);
    ~OffTargetStage() {}
};

class Pipeline {
    std::vector<std::unique_ptr<Stage>> stages_;
public:
    void add(std::unique_ptr<Stage> s) { stages_.push_back(std::move(s)); }
    void run(PipelineContext& ctx);
};