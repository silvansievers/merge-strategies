#ifndef MERGE_AND_SHRINK_MERGE_AND_SHRINK_REPRESENTATION_H
#define MERGE_AND_SHRINK_MERGE_AND_SHRINK_REPRESENTATION_H

#include <cassert>
#include <memory>
#include <vector>

class State;

namespace merge_and_shrink {
class Distances;
class MergeAndShrinkRepresentation {
protected:
    int domain_size;

public:
    explicit MergeAndShrinkRepresentation(int domain_size);
    virtual ~MergeAndShrinkRepresentation() = 0;

    int get_domain_size() const;

    // Store distances instead of abstract state numbers.
    virtual void set_distances(const Distances &) = 0;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) = 0;
    /*
      Return the value that state is mapped to. This is either an abstract
      state (if set_distances has not been called) or a distance (if it has).
      If the represented function is not total, the returned value is DEAD_END
      if the abstract state is PRUNED_STATE or if the (distance) value is INF.
    */
    virtual int get_value(const State &state) const = 0;
    /* Return true iff the represented function is total, i.e., does not map
       to PRUNED_STATE. */
    virtual bool is_total() const = 0;
    virtual bool operator==(const MergeAndShrinkRepresentation &other) const = 0;
    virtual void dump() const = 0;
};


class MergeAndShrinkRepresentationLeaf : public MergeAndShrinkRepresentation {
    const int var_id;

    std::vector<int> lookup_table;
public:
    MergeAndShrinkRepresentationLeaf(int var_id, int domain_size);
    explicit MergeAndShrinkRepresentationLeaf(const MergeAndShrinkRepresentationLeaf *other);
    virtual ~MergeAndShrinkRepresentationLeaf() = default;

    virtual void set_distances(const Distances &) override;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) override;
    virtual int get_value(const State &state) const override;
    virtual bool is_total() const override;
    virtual bool operator==(const MergeAndShrinkRepresentation &other) const override {
        try {
            const MergeAndShrinkRepresentationLeaf &tmp = dynamic_cast<const MergeAndShrinkRepresentationLeaf &>(other);
            assert(domain_size == tmp.domain_size);
            assert(var_id == tmp.var_id);
            assert(lookup_table == tmp.lookup_table);
            return var_id == tmp.var_id && lookup_table == tmp.lookup_table;
        } catch (const std::bad_cast &) {
            assert(false);
            return false;
        }
    }
    virtual void dump() const override;
};


class MergeAndShrinkRepresentationMerge : public MergeAndShrinkRepresentation {
    std::unique_ptr<MergeAndShrinkRepresentation> left_child;
    std::unique_ptr<MergeAndShrinkRepresentation> right_child;
    std::vector<std::vector<int>> lookup_table;
public:
    MergeAndShrinkRepresentationMerge(
        std::unique_ptr<MergeAndShrinkRepresentation> left_child,
        std::unique_ptr<MergeAndShrinkRepresentation> right_child);
    explicit MergeAndShrinkRepresentationMerge(const MergeAndShrinkRepresentationMerge *other);
    virtual ~MergeAndShrinkRepresentationMerge() = default;

    virtual void set_distances(const Distances &distances) override;
    virtual void apply_abstraction_to_lookup_table(
        const std::vector<int> &abstraction_mapping) override;
    virtual int get_value(const State &state) const override;
    virtual bool is_total() const override;
    virtual bool operator==(const MergeAndShrinkRepresentation &other) const override {
        try {
            const MergeAndShrinkRepresentationMerge &tmp = dynamic_cast<const MergeAndShrinkRepresentationMerge &>(other);
            assert(domain_size == tmp.domain_size);
            assert(*left_child.get() == *tmp.left_child.get());
            assert(*right_child.get() == *tmp.right_child.get());
            assert(lookup_table == tmp.lookup_table);
            return *left_child.get() == *tmp.left_child.get() && *right_child.get() == *tmp.right_child.get() && lookup_table == tmp.lookup_table;
        } catch (const std::bad_cast &) {
            assert(false);
            return false;
        }
    }
    virtual void dump() const override;
};
}

#endif
