#include "make_unique.h"
#include "world.h"
#include "quad-tree.h"
#include <algorithm>
#include <iostream>

// TASK 1

// NOTE: You may modify any of the contents of this file, but preserve all function types and names.
// You may add new functions if you believe they will be helpful.



const int QuadTreeLeafSize = 8;
class SequentialNBodySimulator : public INBodySimulator
{
public:
    std::shared_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax)
    {
       // TODO: implement a function that builds and returns a quadtree containing particles.
       // THE PLAN

        //BASE CASE
        // if vector size is 0 throw exception
        if(particles.size()<=0){
            std::cout << "An error occurred, particle vector of size 0 shouldn't be passed into buildQuadTree()" << std::endl;
            return nullptr;
        }
        // if vector size is less than or = to quadtreeleafsize, return a leaf node that holds all the particle in the bmin bmax
        else if (particles.size()<=QuadTreeLeafSize){//this is a leaf node
            auto leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
            leaf->isLeaf=true;//set leaf to true
            leaf->particles= particles;//copy particles
            //(std::static_pointer_cast<QuadTreeNode>(leaf))->isLeaf=true;//set leaf to true
            //(std::static_pointer_cast<QuadTreeNode>(leaf))->particles= particles;//copy particles
            return leaf;
            /*QuadTreeNode leaf;
            leaf.isLeaf = 1;
            leaf.particles = *particles;
            return 
            */
        }

        //RECURSIVE STEP
        // if vector is larger than quadtreeleafsize split vector into 4 vectors for each of the 4 corners of the bmin/bmax bounds
        // and return an empty non leaf node where each child is one of the four vectors
        else{//further divide tree
            //containers for child nodes
            std::vector<Particle> partQuad0;//upper left
            std::vector<Particle> partQuad1;//upper right
            std::vector<Particle> partQuad2;//lower left
            std::vector<Particle> partQuad3;//lower right


            //center
            Vec2 center;
            center.x = (bmin.x+bmax.x);
            center.y = (bmin.y+bmax.y);

            //sort all particles into their particle quadrant
            for(Particle p: particles){
                if(p.position.y > center.y)//quad 0 or 1
                {
                    if(p.position.x < center.x)//left: quad 0
                    {
                        partQuad0.push_back(p);
                    }
                    else//right: quad 1
                    {
                        partQuad1.push_back(p);
                    }
                }
                else//quad 2 or 3
                {
                    if(p.position.x < center.x)//left: quad 2
                    {
                        partQuad2.push_back(p);
                    }
                    else//right: quad 3
                    {
                        partQuad3.push_back(p);
                    }

                }
            }

            //For each Quadrant if particles are greater than 0, create child nodes and add them to the parent node
            auto parent = std::make_shared<QuadTreeNode>();

            //QUAD 0
            if(partQuad0.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = bmin.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = center.x;
                child_bmax.y = bmax.y;

                //create child recursively 
                auto child = buildQuadTree(partQuad0,child_bmin,child_bmax);

                //add child to parent
                parent->children[0] = child;
            }
            
            //QUAD 1
            if(partQuad1.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = bmax.y;

                //create child recursively 
                auto child = buildQuadTree(partQuad1,child_bmin,child_bmax);

                //add child to parent
                parent->children[1] = child;
            }
            
            //QUAD 2
            if(partQuad2.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = bmin.x;
                child_bmin.y = bmin.y;

                Vec2 child_bmax;
                child_bmax.x = center.x;
                child_bmax.y = center.y;

                //create child recursively 
                auto child = buildQuadTree(partQuad2,child_bmin,child_bmax);

                //add child to parent
                parent->children[2] = child;
            }
            
            //QUAD 3
            if(partQuad3.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = bmin.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = center.y;

                
                //create child recursively 
                auto child = buildQuadTree(partQuad3,child_bmin,child_bmax);

                //add child to parent
                
                parent->children[3] = child;
            }
            

            return parent;//a shared pointer to the root node

        }

        return nullptr;
    }
    virtual std::shared_ptr<AccelerationStructure> buildAccelerationStructure(std::vector<Particle> & particles)
    {
        // build quad-tree
        auto quadTree = std::make_shared<QuadTree>();

        // find bounds
        Vec2 bmin(1e30f, 1e30f);
        Vec2 bmax(-1e30f, -1e30f);

        for (auto & p : particles)
        {
            bmin.x = fminf(bmin.x, p.position.x);
            bmin.y = fminf(bmin.y, p.position.y);
            bmax.x = fmaxf(bmax.x, p.position.x);
            bmax.y = fmaxf(bmax.y, p.position.y);
        }

        quadTree->bmin = bmin;
        quadTree->bmax = bmax;

        // build nodes
        quadTree->root = buildQuadTree(particles, bmin, bmax);
        if (!quadTree->checkTree()) {
          std::cout << "Your Tree has Error!" << std::endl;
        }

        return quadTree;
    }
    virtual void simulateStep(AccelerationStructure * accel, std::vector<Particle> & particles, std::vector<Particle> & newParticles, StepParameters params) override
    {
        // TODO: implement sequential version of quad-tree accelerated n-body simulation here,
        // using quadTree as acceleration structure

        //THE PLAN
        //In this case, instead of looking for each particle in a radius, we look for each NODE in a radius
        //We can ignore nodes for which the the closest boundary corner is outside the radius
        //we go through each particle in that nodes and calculate from there
       
    }
};

std::unique_ptr<INBodySimulator> createSequentialNBodySimulator()
{
    return std::make_unique<SequentialNBodySimulator>();
}
