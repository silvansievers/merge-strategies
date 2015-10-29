#ifndef MERGE_AND_SHRINK_HEURISTIC_REPRESENTATION_H
#define MERGE_AND_SHRINK_HEURISTIC_REPRESENTATION_H

#include <cassert>
#include <memory>
#include <vector>

class State;


class HeuristicRepresentation {
protected:
    int domain_size;

public:
    explicit HeuristicRepresentation(int domain_size);
    virtual ~HeuristicRepresentation() = 0;

    int get_domain_size() const;

    virtual int get_abstract_state(const State &state) const = 0;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) = 0;
    virtual bool operator==(const HeuristicRepresentation &other) const = 0;
};


class HeuristicRepresentationLeaf : public HeuristicRepresentation {
    const int var_id;

    std::vector<int> lookup_table;
public:
    HeuristicRepresentationLeaf(int var_id, int domain_size);
    explicit HeuristicRepresentationLeaf(const HeuristicRepresentationLeaf *other);
    virtual ~HeuristicRepresentationLeaf() = default;

    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) override;
    virtual int get_abstract_state(const State &state) const override;
    virtual bool operator==(const HeuristicRepresentation &other) const override {
        try {
            const HeuristicRepresentationLeaf &tmp = dynamic_cast<const HeuristicRepresentationLeaf &>(other);
            assert(domain_size == tmp.domain_size);
            assert(var_id == tmp.var_id);
            assert(lookup_table == tmp.lookup_table);
            return (var_id == tmp.var_id && lookup_table == tmp.lookup_table);
        } catch (const std::bad_cast&) {
            assert(false);
            return false;
        }
    }
};


class HeuristicRepresentationMerge : public HeuristicRepresentation {
    std::unique_ptr<HeuristicRepresentation> left_child;
    std::unique_ptr<HeuristicRepresentation> right_child;
    std::vector<std::vector<int>> lookup_table;
public:
    HeuristicRepresentationMerge(
        std::unique_ptr<HeuristicRepresentation> left_child,
        std::unique_ptr<HeuristicRepresentation> right_child);
    explicit HeuristicRepresentationMerge(const HeuristicRepresentationMerge *other);
    virtual ~HeuristicRepresentationMerge() = default;

    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) override;
    virtual int get_abstract_state(const State &state) const override;
    virtual bool operator==(const HeuristicRepresentation &other) const override {
        try {
            const HeuristicRepresentationMerge &tmp = dynamic_cast<const HeuristicRepresentationMerge &>(other);
            assert(domain_size == tmp.domain_size);
            assert(*left_child.get() == *tmp.left_child.get());
            assert(*right_child.get() == *tmp.right_child.get());
            assert(lookup_table == tmp.lookup_table);
            return (*left_child.get() == *tmp.left_child.get() && *right_child.get() == *tmp.right_child.get() && lookup_table == tmp.lookup_table);
        } catch (const std::bad_cast&) {
            assert(false);
            return false;
        }
    }
};


#endif
