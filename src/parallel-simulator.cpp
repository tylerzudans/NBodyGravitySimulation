#include "make_unique.h"
#include "world.h"
#include "quad-tree.h"
#include <omp.h>
#include <algorithm>
#include <tuple>
#include <iostream>

// TASK 2

// NOTE: You may modify this class definition as you see fit, as long as the class name,
// and type of simulateStep and buildAccelerationStructure remain the same.

const int QuadTreeLeafSize = 8;
class ParallelNBodySimulator : public INBodySimulator
{
public:
    // TODO: implement a function that builds and returns a quadtree containing particles.
    // You do not have to preserve this function type.
    std::shared_ptr<QuadTreeNode> buildQuadTree(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax){
        return buildQuadTreeParallelBasic(particles,bmin,bmax);
    }
    std::shared_ptr<QuadTreeNode> buildQuadTreeParallelBasic(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax){
        //ORCHESTRATION - BUILD TOP OF THE TREE MANUALLY, AND HAVE A SET OF SUBPROBLEMS TO BE DIVIDED AMONGST NODES
        
        //ROOT
        std::shared_ptr<QuadTreeNode> root = std::make_shared<QuadTreeNode>();
        if(particles.size()<=QuadTreeLeafSize) return buildQuadTreeSequential(particles,bmin,bmax);//if tree is small use sequential
        
        //Level 1
        //std::cout << "Dividing into subproblems on level 1 . . ." << "\n";
        std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> subproblems = divideIntoFour(particles,bmin,bmax);
        for(int i = 0; i<subproblems.size(); i++){
            std::shared_ptr<QuadTreeNode> child = buildNode(subproblems[i]);
            if(child->isLeaf) return buildQuadTreeSequential(particles,bmin,bmax);//if tree is small use sequential
            root->children[i]=child;
        }
        //std::cout << "Level 1 children built" << "\n";
        
        #pragma omp parallel num_threads(4) //Build 4 threads and have each take over a child
        {
            int thread_num = omp_get_thread_num();
            auto root_child = root->children[thread_num];
            auto built_child = buildQuadTreeSequential(root_child->particles,std::get<1>(subproblems[thread_num]),std::get<2>(subproblems[thread_num]));
            root->children[omp_get_thread_num()] = built_child;
            //std::cout<<thread_num<<'\n';
        }
        return root;


        /*
        //Level 2
        std::cout << "Dividing into subproblems on level 2 . . ." << "\n";
        subproblems.clear();
        for(auto node: root->children){
            std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> subsubproblems = divideIntoFour(particles,bmin,bmax);
        }
        */


        /*
        while(subproblems.size()<STARTER_NODES){ //build up the subproblems
            std::cout << subproblems.size();
            std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> old_subproblems(subproblems);//make a copy of the old problems
            subproblems.clear();
            for(auto problem: subproblems){//take qeued problems
                for(auto subproblem: divideIntoFour(problem)){//divide each in 4
                    subproblems.push_back(subproblem);//add 4 to queue
                }
                
            }
        }
        */
        std::cout << "Subproblems: " << subproblems.size() << "\n";
        return buildQuadTreeSequential(particles,bmin,bmax);
    }
    std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> divideIntoFour(std::tuple<std::vector<Particle>, Vec2, Vec2> problem_tuple){
        return divideIntoFour(std::get<0>(problem_tuple),std::get<1>(problem_tuple),std::get<2>(problem_tuple));
    }
    std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> divideIntoFour(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax){
        std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>> subdivisions = std::vector<std::tuple<std::vector<Particle>, Vec2, Vec2>>();
        if(particles.size()<=0){
            // leave vector empty
        } 
        else if(particles.size()<=QuadTreeLeafSize){//cannot be subdivided, add function back
            std::tuple<std::vector<Particle>, Vec2, Vec2> quad1 (particles,bmin,bmax);
            subdivisions.push_back(quad1);
        }
        else{//divide into 4 problems
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


            //Create Sub-Problems 
            //if(partQuad0.size()>0){
            if(true){
                Vec2 child_bmin;
                child_bmin = bmin;

                Vec2 child_bmax;
                child_bmax = center;

                std::tuple<std::vector<Particle>, Vec2, Vec2> quad0 (partQuad0,child_bmin,child_bmax);
                subdivisions.push_back(quad0);    
            }
            if(true){
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = bmin.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = center.y;

                std::tuple<std::vector<Particle>, Vec2, Vec2> quad1 (partQuad1,child_bmin,child_bmax);
                subdivisions.push_back(quad1);  
            }
            if(true){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = bmin.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = center.x;
                child_bmax.y = bmax.y;

                std::tuple<std::vector<Particle>, Vec2, Vec2> quad2 (partQuad2,child_bmin,child_bmax);
                subdivisions.push_back(quad2); 
            }
            if(true){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin.x = center.x;
                child_bmin.y = center.y;

                Vec2 child_bmax;
                child_bmax.x = bmax.x;
                child_bmax.y = bmax.y;

                std::tuple<std::vector<Particle>, Vec2, Vec2> quad3 (partQuad3,child_bmin,child_bmax);
                subdivisions.push_back(quad3); 
            }
            
        }
        return subdivisions;

    }
    std::shared_ptr<QuadTreeNode> buildNode(std::tuple<std::vector<Particle>, Vec2, Vec2> problem_tuple){
        return buildNode(std::get<0>(problem_tuple),std::get<1>(problem_tuple),std::get<2>(problem_tuple));
    }
    std::shared_ptr<QuadTreeNode> buildNode(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax){
        if(particles.size()<=0){
            std::cout << "An error occurred, particle vector of size 0 shouldn't be passed into buildQuadTree()" << std::endl;
            return nullptr;
        }
        // if vector size is less than or = to quadtreeleafsize, return a leaf node that holds all the particle in the bmin bmax
        else if (particles.size()<=QuadTreeLeafSize){//this is a leaf node
            std::shared_ptr<QuadTreeNode> leaf = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
            leaf->isLeaf=true;//set leaf to true
            leaf->particles= particles;//copy particles
            return leaf;
        }
        else{
            std::shared_ptr<QuadTreeNode> node = std::make_shared<QuadTreeNode>();//leaf is of type shared_ptr<QuadTreeNode>
            //leaf->isLeaf=true;//set leaf to true
            node->particles= particles;//copy particles
            return node;
        }
    }
    std::shared_ptr<QuadTreeNode> buildQuadTreeSequential(std::vector<Particle> & particles, Vec2 bmin, Vec2 bmax)
    {
        //BASE CASE
        // if vector size is 0 throw exception
        if(particles.size()<=0){
            std::cout << "An error occurred, particle vector of size 0 shouldn't be passed into buildQuadTree()" << std::endl;
            return nullptr;
        }
        // if vector size is less than or = to quadtreeleafsize, return a leaf node that holds all the particle in the bmin bmax
        else if (particles.size()<=QuadTreeLeafSize){//this is a leaf node
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

            //std::cout << "Dividing particles into 4 quadrants . . ." << "\n";

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
            std::shared_ptr<QuadTreeNode> parent = std::make_shared<QuadTreeNode>();

            //std::cout << "Recursively building children . . ." << "\n";

            //QUAD 0 - nw
            if(partQuad0.size()>0){//child must be created in quad 0
                //set child boundaries
                Vec2 child_bmin;
                child_bmin = bmin;

                Vec2 child_bmax;
                child_bmax = center;

                //create child recursively 
                std::shared_ptr<QuadTreeNode> child = buildQuadTreeSequential(partQuad0,child_bmin,child_bmax);

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
                std::shared_ptr<QuadTreeNode> child = buildQuadTreeSequential(partQuad1,child_bmin,child_bmax);

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
                std::shared_ptr<QuadTreeNode> child = buildQuadTreeSequential(partQuad2,child_bmin,child_bmax);

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
                std::shared_ptr<QuadTreeNode> child = buildQuadTreeSequential(partQuad3,child_bmin,child_bmax);

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

    // Do not modify this function type.
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

    // Do not modify this function type.
    virtual void simulateStep(AccelerationStructure * accel,
                            std::vector<Particle> & particles,
                            std::vector<Particle> & newParticles,
                            StepParameters params) override
    {
        
       #pragma omp parallel for
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

// Do not modify this function type.
std::unique_ptr<INBodySimulator> createParallelNBodySimulator()
{
  return std::make_unique<ParallelNBodySimulator>();
}
