#include "bbox.h"

class Geometry;
class IntersectionInfo;

class BVHTree {
private:
	struct PrimInfo
	{
		PrimInfo(size_t idx, BBox bounds) : primitiveIdx(idx), boundingBox(bounds), centroid(.5f * bounds.vmin + .5f * bounds.vmax)
		{
		}
		size_t primitiveIdx;
		BBox boundingBox;
		Vector centroid;
	};

	struct Node
	{
		void initLeaf(int first, int n, const BBox& b)
		{
			firstPrimOffset = first;
			primitiveCount = n;
			bounds = b;
			children[0] = children[1] = nullptr;
		}

		void initInterior(int axis, Node* child1, Node* child2)
		{
			bounds.extend(child1->bounds);
			bounds.extend(child2->bounds);
			splitAxis = axis;
			primitiveCount = 0;
			children[0] = child1;
			children[1] = child2;
		}

		int primitiveCount, firstPrimOffset, splitAxis;
		BBox bounds;
		Node* children[2];
	};

	struct Treelet
	{
		int startIdx, primitiveCount;
		Node* nodes;
	};

	struct MortonPrim {
		int primitiveIndex;
		uint64_t mortonCode;
	};

	struct LinearNode
	{
		BBox bounds;
		union
		{
			int primitivesOffset;
			int secondChildOffset;
		};

		uint16_t primitiveCount;
		uint8_t axis; // axis interior nodes were split on
		// uint8_t pad[1]; // padding for 32b... big bounding boxes... 56 bytes
	};
public:
	~BVHTree();
	void addGeometry(Geometry *geometry);
	void build();
	bool intersect(const Ray& ray, float tMin, float tMax, IntersectionInfo& intersection);
private:
	void clear();
	uint64_t weirdShift(uint64_t x);
	uint64_t encodeMorton3(const Vector& val);
	Node* buildTreelets(Node *&buildNodes, MortonPrim* mortonPrims, int primitiveCount, int& totalNodes, int& orderedPrimsOffset, int bitIdx);
	Node* connectTreelets(std::vector<Node*>& roots, int start, int end, int& totalNodes) const;
	int flatten(Node* node, int& offset);

	std::vector<PrimInfo> m_Primitives;
	std::vector<Geometry*> m_OrderedPrims;
	std::vector<Geometry*> m_FinalPrims;
	LinearNode* m_SearchNodes = nullptr;
	uint32_t m_MaxPrimsPerNode = 1;
	float m_IntersectionCost = 1.0f; // cost of calculating intersection

	uint32_t m_PrimIdx = 0;
};
