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

        //BASE CASE
        // if vector size is 0 throw exception
        if(particles.size()<=0){
            std::cout << "An error occurred, particle vector of size 0 shouldn't be passed into buildQuadTree()" << std::endl;
            return nullptr;
        }
        // if vector size is less than or = to quadtreeleafsize, return a leaf node that holds all the particle in the bmin bmax
        else if (particles.size()<=QuadTreeLeafSize){//this is a leaf node
            /*
            //DEBUG 
            std::cout<<"Quad with particle vector of size "<<particles.size()<<" -> ("<<bmin.x<<','<<bmin.y<<") - ("<<bmax.x<<','<<bmax.y<<")"<<std::endl; 
            for(Particle p : particles){
                std::cout<<p.id<<" -> "<<p.position.x<<','<<p.position.y<<std::endl;
            }
            std::cout<<std::endl;
            */

            std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
            leaf->isLeaf=true;//set leaf to true
            leaf->particles= particles;//copy particles
            return leaf;
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
            center.x = (bmin.x+bmax.x)/2;
            center.y = (bmin.y+bmax.y)/2;

            //sort all particles into their particle quadrant
            for(Particle p: particles){
                if(p.position.y <  center.y)//quad 0 or 1
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
            //auto parent = std::make_shared<QuadTreeNode>();
            std::shared_ptr<QuadTreeNode> parent = std::make_shared<QuadTreeNode>();

            //QUAD 0 - nw
            if(partQuad0.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin = bmin;

                Vec2 child_bmax;
                child_bmax = center;

                //create child recursively 
                std::shared_ptr<QuadTreeNode> child = buildQuadTree(partQuad0,child_bmin,child_bmax);

                //add child to parent
                parent->children[0] = child;
            }
            else{//child is empty leaf node
                std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
                leaf->isLeaf=true;//set leaf to true
                parent->children[0] = leaf;
            }
            
            //QUAD 1 - ne
            if(partQuad1.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = bmin.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = center.y;

                //create child recursively 
                std::shared_ptr<QuadTreeNode> child = buildQuadTree(partQuad1,child_bmin,child_bmax);

                //add child to parent
                parent->children[1] = child;
            }
            else{//child is empty leaf node
                std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
                leaf->isLeaf=true;//set leaf to true
                parent->children[1] = leaf;
            }
            
            //QUAD 2 -sw
            if(partQuad2.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = bmin.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = center.x;
                child_bmax.y = bmax.y;

                //create child recursively 
                std::shared_ptr<QuadTreeNode> child = buildQuadTree(partQuad2,child_bmin,child_bmax);

                //add child to parent
                parent->children[2] = child;
            }
            else{//child is empty leaf node
                std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
                leaf->isLeaf=true;//set leaf to true
                parent->children[2] = leaf;
            }
            
            //QUAD 3 -se
            if(partQuad3.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = bmax.y;

                
                //create child recursively 
                std::shared_ptr<QuadTreeNode> child = buildQuadTree(partQuad3,child_bmin,child_bmax);

                //add child to parent
                
                parent->children[3] = child;
            }
            else{//child is empty leaf node
                std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
                leaf->isLeaf=true;//set leaf to true
                parent->children[3] = leaf;
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
       //#pragma omp parallel for
        for (int i = 0; i < (int)particles.size(); i++)
        {
            auto pi = particles[i];
            Vec2 force = Vec2(0.0f, 0.0f);
            std::vector<Particle> particlesInRange;
            accel->getParticles(particlesInRange,pi.position,params.cullRadius);//get proximal oarticles
            for(Particle pr : particlesInRange) force += computeForce(pi,pr,params.cullRadius);//add force to particle


            // update particle state using the computed force
            newParticles[i] = updateParticle(pi, force, params.deltaTime);
        }
    }
};

std::unique_ptr<INBodySimulator> createSequentialNBodySimulator()
{
    return std::make_unique<SequentialNBodySimulator>();
}
