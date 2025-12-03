#ifndef NODE_H
#define NODE_H

#include <iostream>

/**
 * Node structure representing an element in the Disjoint Set.
 * Templated to accept any type T that supports comparison (<, ==) and stream output.
 */
template <typename T>
struct Node {
    T value;
    Node* parent; // nullptr implies this node is a root

    Node(T v) : value(v), parent(nullptr) {}
};

/**
 * Overload << operator for easy printing of Nodes.
 * Note: Must be inline or template to be in header.
 */
template <typename T>
std::ostream& operator<<(std::ostream& os, const Node<T>& node) {
    os << "Node(value=" << node.value << ", parent=";
    if (node.parent) {
        os << node.parent->value;
    } else {
        os << "None";
    }
    os << ")";
    return os;
}

#endif // NODE_H