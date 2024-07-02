#include "bbox.h"
#include "bvh.h"
#include "geometry.h"

#include <functional>

BVHTree::~BVHTree() {
	clear();
}

void BVHTree::addGeometry(Geometry *geometry)
{
	BBox box;
	geometry->expandBox(box);
	m_Primitives.push_back({ m_PrimIdx++, box });
	m_FinalPrims.push_back(geometry);
}

void BVHTree::clear()
{
	delete[] m_SearchNodes;
	m_SearchNodes = nullptr;
}

uint64_t BVHTree::weirdShift(uint64_t x) // pbr book, but this is with 64 bits
{
	x = (x | (x << 32)) & 0x001f00000000ffff; // 0000000000011111000000000000000000000000000000001111111111111111
	x = (x | (x << 16)) & 0x001f0000ff0000ff; // 0000000000011111000000000000000011111111000000000000000011111111
	x = (x | (x <<  8)) & 0x100f00f00f00f00f; // 0001000000001111000000001111000000001111000000001111000000000000
	x = (x | (x <<  4)) & 0x10c30c30c30c30c3; // 0001000011000011000011000011000011000011000011000011000100000000
	x = (x | (x <<  2)) & 0x1249249249249249; // 0001001001001001001001001001001001001001001001001001001001001001
	return x;
}

uint64_t BVHTree::encodeMorton3(const Vector& val)
{
	return (weirdShift(val.z) << 2) | (weirdShift(val.y) << 1) | weirdShift(val.x);
}

void BVHTree::build()
{
	{
		m_MaxPrimsPerNode = 1;
		m_IntersectionCost = 2.0f;
	}
	// Timer timer;
	printf("Building BVH with %d primitives\n", (int)m_Primitives.size());
	BBox bounds;
	for (const auto& prim : m_Primitives) // Bounding box of all primitives
		bounds.add(prim.centroid);

	std::vector<MortonPrim> mortonPrims;
	mortonPrims.resize(m_Primitives.size());
	const int mortonBits = 21; // so we can use 21 bits for each axis with int = 3x21 63
	const int mortonScale = 1 << mortonBits; // Multiply by 2^21 since I can fit 21 bits in the morton thing
	for (int i = 0; i < m_Primitives.size(); i++) // pbr book does this in parallel
	{
		mortonPrims[i].primitiveIndex = m_Primitives[i].primitiveIdx;
		Vector centroidOffset = bounds.offset(m_Primitives[i].centroid);
		mortonPrims[i].mortonCode = encodeMorton3(centroidOffset * mortonScale);
	}

	std::sort(mortonPrims.begin(), mortonPrims.end(), [](const MortonPrim& l, const MortonPrim& r){ return l.mortonCode < r.mortonCode; }); // pbr book uses radix sort here

	std::vector<Treelet> treeletsToBuild;
	int start = 0;
	for (int end = 1; end < (int)mortonPrims.size(); end++)
	{
		uint64_t mask = 0x3ffc000000000000; // top 12 bits, divide them into groups whose top 12 bits match
		if ((mortonPrims[start].mortonCode & mask) != (mortonPrims[end].mortonCode & mask))
		{
			int primitiveCount = end - start;
			int maxBVHNodes = 2 * primitiveCount;
			Node* nodes = new Node[maxBVHNodes];
			treeletsToBuild.push_back({ start, primitiveCount, nodes });
			start = end;
		}
	}

	printf("%lld treelets", treeletsToBuild.size() + 1);

	int primitiveCount = mortonPrims.size() - start;
	int maxBVHNodes = 2 * primitiveCount;
	treeletsToBuild.push_back({ start, primitiveCount, new Node[maxBVHNodes] });

	// Could also do this in parallel
	int orderedPrimsOffset = 0;
	int totalNodes = 0;
	const int firstBitIndex = 62 - 12;
	m_OrderedPrims.resize(m_Primitives.size());
	for (int i = 0; i < treeletsToBuild.size(); i++)
		treeletsToBuild[i].nodes = buildTreelets(treeletsToBuild[i].nodes, &mortonPrims[treeletsToBuild[i].startIdx], treeletsToBuild[i].primitiveCount, totalNodes, orderedPrimsOffset, firstBitIndex);

	std::vector<Node*> finishedTreelets; // Create the rest of the tree using SAH
	finishedTreelets.reserve(treeletsToBuild.size());
	for (Treelet& treelet : treeletsToBuild)
		finishedTreelets.push_back(treelet.nodes);
	Node* root = connectTreelets(finishedTreelets, 0, finishedTreelets.size(), totalNodes);
	m_SearchNodes = new LinearNode[totalNodes];
	m_FinalPrims.swap(m_OrderedPrims);
	m_Primitives.clear();

	std::function<void(Node*, const std::string&, bool)> printBVH = [&](Node* node, const std::string& prefix, bool isLeft) {
		printf("%s", prefix.c_str());
		printf(isLeft ? "|--" : "L--");
		if (node->children[0] != nullptr)
			printf("Interior: %f, %f, %f, %f, %f, %f\n", node->bounds.vmin.x, node->bounds.vmin.y, node->bounds.vmin.z, node->bounds.vmax.x, node->bounds.vmax.y, node->bounds.vmax.z);
		else
			printf("Leaf: %f, %f, %f, %f, %f, %f\n", node->bounds.vmin.x, node->bounds.vmin.y, node->bounds.vmin.z, node->bounds.vmax.x, node->bounds.vmax.y, node->bounds.vmax.z);
		if (node->children[0] != nullptr)
			printBVH(node->children[0], prefix + (isLeft ? "|   " : "    "), true);
		if (node->children[1] != nullptr)
			printBVH(node->children[1], prefix + (isLeft ? "|   " : "    "), false);
	};
	// printBVH(root, "", false);
	int32_t offset = 0;
	flatten(root, offset); // pbr book
	// LOG_ACCEL_BUILD(AcceleratorType::BVH, timer.toMs<float>(timer.elapsedNs() / 1000.0f), totalNodes, totalNodes * sizeof(LinearNode) + sizeof(*this) + sizeof(m_FinalPrims[0]) * m_FinalPrims.size());
	// printf("Built BVH with %d nodes in %f seconds\n", totalNodes, Timer::toMs<float>(timer.elapsedNs()) / 1000.0f);
}

BVHTree::Node* BVHTree::buildTreelets(Node *&buildNodes, MortonPrim* mortonPrims, int primitiveCount, int& totalNodes, int& orderedPrimsOffset, int bitIdx)
{
	if (bitIdx == -1 || primitiveCount < m_MaxPrimsPerNode) // We need to create a leaf, either because we can fit the nodes left in a single leaf, or because we can't split
	{
		totalNodes++;
		Node* node = buildNodes++;
		BBox bounds;
		int firstPrimOffset = orderedPrimsOffset;
		orderedPrimsOffset += primitiveCount;
		for (int i = 0; i < primitiveCount; i++)
		{
			int primitiveIdx = mortonPrims[i].primitiveIndex;
			m_OrderedPrims[firstPrimOffset + i] = m_FinalPrims[primitiveIdx];
			bounds.extend(m_Primitives[primitiveIdx].boundingBox);
		}
		node->initLeaf(firstPrimOffset, primitiveCount, bounds);
		return node;
	}
	else // Create an internal node with two children
	{
		int mask = 1 << bitIdx;
		if ((mortonPrims[0].mortonCode & mask) == ((mortonPrims[primitiveCount - 1].mortonCode & mask))) // Check if all nodes are on the same side of the plane
			return buildTreelets(buildNodes, mortonPrims, primitiveCount, totalNodes, orderedPrimsOffset, bitIdx - 1);
		int l = 0, r = primitiveCount - 1;
		while (l + 1 != r) // binary search for region where bitIdx bit goes from 0 to 1
		{
			int mid = (l + r) / 2;
			// std::cout << std::bitset<32>(mortonPrims[l].mortonCode).to_string() << std::endl;
			// std::cout << std::bitset<32>(mortonPrims[mid].mortonCode).to_string() << std::endl;
			if ((mortonPrims[l].mortonCode & mask) == (mortonPrims[mid].mortonCode & mask))
				l = mid;
			else
				r = mid;
		}
		int splitOffset = r; // Primitives are already on correct sides of plane
		totalNodes++;
		Node* node = buildNodes++;
		Node* node1 = buildTreelets(buildNodes, mortonPrims, splitOffset, totalNodes, orderedPrimsOffset, bitIdx - 1);
		Node* node2 = buildTreelets(buildNodes, &mortonPrims[splitOffset], primitiveCount - splitOffset, totalNodes, orderedPrimsOffset, bitIdx - 1);
		int axis = bitIdx % 3;
		node->initInterior(axis, node1, node2);
		return node;
	}
}

int BVHTree::flatten(Node* node, int& offset)
{
	// Store the tree in dfs parent left right order
	LinearNode* linearNode = &m_SearchNodes[offset];
	linearNode->bounds = node->bounds;
	int myOffset = offset++;
	if (node->primitiveCount > 0)
	{
		linearNode->primitivesOffset = node->firstPrimOffset;
		linearNode->primitiveCount = node->primitiveCount;
	}
	else
	{
		linearNode->axis = node->splitAxis;
		linearNode->primitiveCount = 0;
		flatten(node->children[0], offset);
		linearNode->secondChildOffset = flatten(node->children[1], offset);
	}
	return myOffset;
}

bool BVHTree::intersect(const Ray& ray, float tMin, float tMax, IntersectionInfo& intersection)
{
	Vector invDir = ray.dir.reciprocal();
	int negativeDir[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
	// Offset of next element in stack, offset in nodes list
	int toVisitOffset = 0, currentNodeIndex = 0;
	int nodesToVisit[64];
	bool hit = false;
	while (true)
	{
		const LinearNode* node = &m_SearchNodes[currentNodeIndex];
		if (node->bounds.testIntersect(ray))
		{
			if (node->primitiveCount > 0) // leaf
			{
				for (int i = 0; i < node->primitiveCount; i++)
				{
					if (m_FinalPrims[node->primitivesOffset + i]->intersect(ray, intersection, tMin, tMax))
					{
						// return true;
						hit = true; // Need to keep going, since there might be closer intersections, so just update
						tMax = intersection.dist;
					}
				}

				if (toVisitOffset == 0) break;
				currentNodeIndex = nodesToVisit[--toVisitOffset];
			}
			else // interior, so visit child
			{
				if (negativeDir[node->axis]) // if the axis we split on has negative direction, visit the second child. In 2D this is:
				{
					/*
					Let's say we split on the x axis. Then we want to visit the second child first if the ray is going from right to left, and the first child if the ray is going from left to right.
					This way we can easily discard the second child, since it would hit the box on the left and it's closer

			    *
				 \
				  \
				   \
				    >

					---------   |
					|       |   |
					|       |   |
					---------   |     ---------
								|     |       |
								|     |       |
								|     ---------
								|
					*/
					nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
					currentNodeIndex = node->secondChildOffset;
				}
				else
				{
					nodesToVisit[toVisitOffset++] = node->secondChildOffset;
					currentNodeIndex++;
				}
			}
		}
		else
		{
			if (toVisitOffset == 0)
				break;
			currentNodeIndex = nodesToVisit[--toVisitOffset];
		}
	}
	return hit;
}

BVHTree::Node* BVHTree::connectTreelets(std::vector<Node*>& roots, int start, int end, int& totalNodes) const {
	int nodeCount = end - start;
	if (nodeCount== 1) return roots[start];
	totalNodes++;
	Node* node = new Node();
	BBox bounds;
	for (int i = start; i < end; i++)
		bounds.extend(roots[i]->bounds);

	// Only considering the centroids of objects. With these scenes the objects are pretty much the same size, but with differently sized objects this wouldn't work as well
	BBox centroidBounds;
	for (int i = start; i < end; i++)
	{
		Vector centroid = (roots[i]->bounds.vmin + roots[i]->bounds.vmax) * 0.5f;
		centroidBounds.add(centroid);
	}

	int dim = centroidBounds.maxExtent(); // Divide on largest axis, Maybe worth checking all 3?
	const int bucketCount = 12; // Put everything in buckets and try to cut between the buckets. Choose the one with the best cost
	struct BucketInfo {
		int count = 0;
		BBox bounds;
	};
	BucketInfo buckets[16];
	for (int i = start; i < end; i++) // Put into buckets
	{
		float centroid = (roots[i]->bounds.vmin[dim] + roots[i]->bounds.vmax[dim]) * 0.5f;
		int b = bucketCount * (((centroid - centroidBounds.vmin[dim]) / (centroidBounds.vmax[dim] - centroidBounds.vmin[dim])));
		if (b >= bucketCount)
			b = bucketCount - 1;
		buckets[b].count++;
		buckets[b].bounds.extend(roots[i]->bounds);
		// printf("%f, %d %d\n", centroid, b);
	}

	const float traversalCost = 0.125f; // cost of figuring out which child to visit

	float cost[bucketCount - 1];
	// printf("Start: %d, end: %d, dim: %d\n", start, end, dim);
	// printf("Centroid: %f, %f, %f\n", centroidBounds.min.x, centroidBounds.min.y, centroidBounds.min.z);
	// printf("Centroid: %f, %f, %f\n", centroidBounds.max.x, centroidBounds.max.y, centroidBounds.max.z);
	for (int i = 0; i < bucketCount - 1; i++)
	{
		// printf("Bucket %d: %d\n", i, buckets[i].count);
		BBox b0, b1;
		int count0 = 0, count1 = 0;
		for (int j = 0; j <= i; j++)
		{
			b0.extend(buckets[j].bounds);
			count0 += buckets[j].count;
		}
		for (int j = 0; j <= i; j++)
		{
			b1.extend(buckets[j].bounds);
			count1 += buckets[j].count;
		}
		cost[i] = traversalCost + m_IntersectionCost * (count0 * b0.area() + count1 * b1.area()) / bounds.area();
	}

	float minCost = cost[0];
	int minCostBucketIdx = 0;
	for (int i = 1; i < bucketCount - 1; i++)
	{
		if (cost[i] < minCost)
		{
			minCost = cost[i];
			minCostBucketIdx = i;
		}
	}

	// 100% not copied from pbr book
	Node** pmid = std::partition(&roots[start], &roots[end - 1] + 1, [=](const Node* node)
		{
			float centroid = (node->bounds.vmin[dim] + node->bounds.vmax[dim]) * 0.5f;
			int b = bucketCount * (((centroid - centroidBounds.vmin[dim]) / (centroidBounds.vmax[dim] - centroidBounds.vmin[dim])));
			if (b >= bucketCount)
				b = bucketCount - 1;
			return b <= minCostBucketIdx;
		});
	int mid = pmid - &roots[0];
	node->initInterior(dim, connectTreelets(roots, start, mid, totalNodes), connectTreelets(roots, mid, end, totalNodes));
	return node;
}