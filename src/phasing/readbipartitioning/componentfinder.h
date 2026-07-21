#ifndef COMPONENT_FINDER_H
#define COMPONENT_FINDER_H

#include <vector>
#include <map>
#include <iostream>
#include "node.h"

/**
 * ComponentFinder class
 * Implements Union-Find where the representative is always the smallest value in the set.
 */
template <typename T>
class ComponentFinder {
    private:
        // We use std::map because it guarantees pointer stability.
        std::map<T, Node<T>> nodes;

        /**
         * Internal helper to find the root node with path compression.
         */
        Node<T>* _find_node(const T& value) {
            Node<T>* node = &nodes.at(value);
            Node<T>* root = node;

            // Find root
            while (root->parent != nullptr) {
                root = root->parent;
            }

            // Path compression
            while (node->parent != nullptr) {
                Node<T>* next = node->parent;
                node->parent = root;
                node = next;
            }

            return root;
        }

public:
    /**
     * Constructor: Initializes the sets with the provided values.
     */
    ComponentFinder(const std::vector<T>& values) {
        for (const auto& val : values) {
            nodes.emplace(val, Node<T>(val));
        }
    }

    /**
     * Merges the sets containing x and y.
     */
    void merge(const T& x, const T& y) {
        if (x == y) return;

        Node<T>* x_root = _find_node(x);
        Node<T>* y_root = _find_node(y);

        if (x_root == y_root) {
            return;
        }

        // Merge smaller value as parent
        if (x_root->value < y_root->value) {
            y_root->parent = x_root;
        } else {
            x_root->parent = y_root;
        }
    }

    /**
     * Return which component x belongs to.
     */
    T find(const T& value) {
        return _find_node(value)->value;
    }

    /**
     * Returns the total number of connected components.
     */
    size_t size() {
        size_t count = 0;
        for (const auto& pair : nodes) {
            // A node is a component representative (root) if its parent is null
            if (pair.second.parent == nullptr) {
                ++count;
            }
        }
        return count;
    }

    /**
     * Debug print function.
     */
    void print() {
        for (auto& pair : nodes) {
            T val = pair.first;
            Node<T>& node = pair.second;
            Node<T>* root = _find_node(val);
            
            std::cout << val << " : " << node 
                      << " is represented by " << root->value << std::endl;
        }
    }
};

#endif