/*
 * Buddy Page Allocation Algorithm
 * SKELETON IMPLEMENTATION -- TO BE FILLED IN FOR TASK (2)
 */

/*
 * STUDENT NUMBER: s1758009
 */
#include <infos/mm/page-allocator.h>
#include <infos/mm/mm.h>
#include <infos/kernel/kernel.h>
#include <infos/kernel/log.h>
#include <infos/util/math.h>
#include <infos/util/printf.h>

using namespace infos::kernel;
using namespace infos::mm;
using namespace infos::util;

#define MAX_ORDER	17

/**
 * A buddy page allocation algorithm.
 */
class BuddyPageAllocator : public PageAllocatorAlgorithm
{
private:
	/**
	 * Returns the number of pages that comprise a 'block', in a given order.
	 * @param order The order to base the calculation off of.
	 * @return Returns the number of pages in a block, in the order.
	 */
	static inline constexpr uint64_t pages_per_block(int order)
	{
		/* The number of pages per block in a given order is simply 1, shifted origin by the order number.
		 * For example, in order-2, there are (1 << 2) == 4 pages in each block.
		 */
		return (1 << order);
	}
	
	/**
	 * Returns TRUE if the supplied page descriptor is correctly aligned for the 
	 * given order.  Returns FALSE otherwise.
	 * @param pgd The page descriptor to test alignment for.
	 * @param order The order to use for calculations.
	 */
	static inline bool is_correct_alignment_for_order(const PageDescriptor *pgd, int order)
	{
		// Calculate the page-frame-number for the page descriptor, and return TRUE if
		// it divides evenly into the number pages in a block of the given order.
		return (sys.mm().pgalloc().pgd_to_pfn(pgd) % pages_per_block(order)) == 0;
	}
	
	/** Given a page descriptor, and an order, returns the buddy PGD.  The buddy could either be
	 * to the origin or the buddy of PGD, in the given order.
	 * @param pgd The page descriptor to find the buddy for.
	 * @param order The order in which the page descriptor lives.
	 * @return Returns the buddy of the given page descriptor, in the given order.
	 */
	PageDescriptor *buddy_of(PageDescriptor *pgd, int order)
	{
		// (1) Make sure 'order' is within range
		if (order >= MAX_ORDER) {
			return NULL;
		}

		// (2) Check to make sure that PGD is correctly aligned in the order
		if (!is_correct_alignment_for_order(pgd, order)) {
			return NULL;
		}
				
		// (3) Calculate the page-frame-number of the buddy of this page.
		// * If the PFN is aligned to the next order, then the buddy is the next block in THIS order.
		// * If it's not aligned, then the buddy must be the previous block in THIS order.
		uint64_t buddy_pfn = is_correct_alignment_for_order(pgd, order + 1) ?
			sys.mm().pgalloc().pgd_to_pfn(pgd) + pages_per_block(order) : 
			sys.mm().pgalloc().pgd_to_pfn(pgd) - pages_per_block(order);
		
		// (4) Return the page descriptor associated with the buddy page-frame-number.
		return sys.mm().pgalloc().pfn_to_pgd(buddy_pfn);
	}
	
	/**
	 * Inserts a block into the free list of the given order.  The block is inserted in ascending order.
	 * @param pgd The page descriptor of the block to insert.
	 * @param order The order in which to insert the block.
	 * @return Returns the slot (i.e. a pointer to the pointer that points to the block) that the block
	 * was inserted into.
	 */
	PageDescriptor **insert_block(PageDescriptor *pgd, int order)
	{
		// Starting from the _free_area array, find the slot in which the page descriptor
		// should be inserted.
		PageDescriptor **slot = &_free_areas[order];
		
		// Iterate whilst there is a slot, and whilst the page descriptor pointer is numerically
		// greater than what the slot is pointing to.
		while (*slot && pgd > *slot) {
			slot = &(*slot)->next_free;
		}
		
		// Insert the page descriptor into the linked list.
		pgd->next_free = *slot;
		*slot = pgd;
		
		// Return the insert point (i.e. slot)
		return slot;
	}
	
	/**
	 * Removes a block from the free list of the given order.  The block MUST be present in the free-list, otherwise
	 * the system will panic.
	 * @param pgd The page descriptor of the block to remove.
	 * @param order The order in which to remove the block from.
	 */
	void remove_block(PageDescriptor *pgd, int order)
	{
		// Starting from the _free_area array, iterate until the block has been located in the linked-list.
		PageDescriptor **slot = &_free_areas[order];
		while (*slot && pgd != *slot) {
			slot = &(*slot)->next_free;
		}

		

		// Make sure the block actually exists.  Panic the system if it does not.
		assert(*slot == pgd);
		
		// Remove the block from the free list.
		*slot = pgd->next_free;
		pgd->next_free = NULL;
	}
	
	/**
	 * Given a pointer to a block of free memory in the order "source_order", this function will
	 * split the block in half, and insert it into the order below.
	 * @param block_pointer A pointer to a pointer containing the beginning of a block of free memory.
	 * @param source_order The order in which the block of free memory exists.  Naturally,
	 * the split will insert the two new blocks into the order below.
	 * @return Returns the origin-hand-side of the new block.
	 */
	PageDescriptor *split_block(PageDescriptor **block_pointer, int source_order)
	{
		// Make sure there is an incoming pointer.
		assert(*block_pointer);
		
		// Make sure the block_pointer is correctly aligned.
		assert(is_correct_alignment_for_order(*block_pointer, source_order));

		// Make sure the source_order is in suitable range.
		assert(source_order >= 0 && (source_order <= (MAX_ORDER-1)));

		// Can't split down if source order is zero, so directly return.
		if (source_order == 0){
			return *block_pointer; 
		}
		
		// Set the target block size and order that you want to split to.	
		int target_order = source_order - 1;
		int target_block_size = pages_per_block(target_order);
		
		// Since split happens on source order, the origin page descriptor be the new left block and new right block must be its buddy.
		// Thus we could directly calculate the offset to get the new right block.
		PageDescriptor* new_left = *block_pointer;
		PageDescriptor* new_right = new_left + target_block_size;

		// Make sure the products are buddy.
		assert(buddy_of(new_left, target_order) == new_right);

		// Make sure the products are correctly aligned.
		assert(is_correct_alignment_for_order(new_right, target_order));
		assert(is_correct_alignment_for_order(new_left, target_order));

		// Remove the old block from _free_areas and add new blocks.
		remove_block(new_left, source_order);
		insert_block(new_left, target_order);
		insert_block(new_right, target_order);

		return new_left;
	}
	
	/**
	 * Takes a block in the given source order, and merges it (and it's buddy) into the next order.
	 * This function assumes both the source block and the buddy block are in the free list for the
	 * source order.  If they aren't this function will panic the system.
	 * @param block_pointer A pointer to a pointer containing a block in the pair to merge.
	 * @param source_order The order in which the pair of blocks live.
	 * @return Returns the new slot that points to the merged block.
	 */
	PageDescriptor **merge_block(PageDescriptor **block_pointer, int source_order)
	{
		// Make sure the argument is valid.
		assert(*block_pointer);
		
		// Make sure the area_pointer is correctly aligned.
		assert(is_correct_alignment_for_order(*block_pointer, source_order));

		// Make sure the source_order is in suitable range.
		assert(source_order >= 0 && (source_order <= (MAX_ORDER-1)));

		// Can't merge up if source order is up to maximum, so directly return.
		if (source_order == MAX_ORDER-1){
			return block_pointer; 
		}

		// Find the buddy block that you want to merge and target order.
		PageDescriptor* to_merge = buddy_of(*block_pointer, source_order);
		int target_order = source_order + 1;

		// Remove the old blocks.
		remove_block(*block_pointer, source_order);
		remove_block(to_merge, source_order);

		// Declare the pgd of new block.
		PageDescriptor* product = NULL;

		// The argument block could be left of right in a pair of buddy
		// Since they "sum" to a new block, the page descriptor of left block in buddy should always be pgd of new block.
		// Decide which block is left block and assign the new block.
		if (*block_pointer > to_merge) {
			product = to_merge;
		}else{
			product = *block_pointer;
		}

		// Make sure the new block is correctly aligned.
		assert(is_correct_alignment_for_order(product, target_order));

		// Insert.
		PageDescriptor** to_return = insert_block(product, target_order);
		
		return to_return;		
	}
	
public:
	/**
	 * Constructs a new instance of the Buddy Page Allocator.
	 */
	BuddyPageAllocator() {
		// Iterate over each free area, and clear it.
		for (unsigned int i = 0; i < ARRAY_SIZE(_free_areas); i++) {
			_free_areas[i] = NULL;
		}
	}
	
	/**
	 * Allocates 2^order number of contiguous pages
	 * @param order The power of two, of the number of contiguous pages to allocate.
	 * @return Returns a pointer to the first page descriptor for the newly allocated page range, or NULL if
	 * allocation failed.
	 */
	PageDescriptor *alloc_pages(int order) override
	{
		// Make sure the source_order is in suitable range.
		assert(order >= 0 && (order <= (MAX_ORDER-1)));

		// If there is available slot in this order, directly return it.
		if (_free_areas[order]) {
			PageDescriptor* to_return = _free_areas[order];
			remove_block(to_return, order);
			return to_return;
		}else{
			// Otherwise, check the _free_areas up to find available slot.
			for(int i = order; i < MAX_ORDER; i++)
			{
				// When find an available slot: need to split it down to requested order and allocate.
				if (_free_areas[i]) {
					PageDescriptor* to_split = _free_areas[i];

					// Iteratively split down to requested order, and only return the "first" block in these split operation.
					while(i > order){
						to_split = split_block(&to_split, i);
						i -= 1;
					}

					// Make sure the product block is correctly aligned to requested order.
					assert(is_correct_alignment_for_order(to_split, order));

					remove_block(to_split, order);
					return to_split;
				}
			}
		}

		// If function executes to here, it means the function can't find am available block. It fails and returns NULL. 
		return NULL;
	}
	
	/**
	 * Frees 2^order contiguous pages.
	 * @param pgd A pointer to an array of page descriptors to be freed.
	 * @param order The power of two number of contiguous pages to free.
	 */
	void free_pages(PageDescriptor *pgd, int order) override
	{
		// Make sure that the incoming page descriptor is correctly aligned
		// for the order on which it is being freed, for example, it is
		// illegal to free page 1 in order-1.
		assert(is_correct_alignment_for_order(pgd, order));

		// Insert the block back to _free_areas.
		insert_block(pgd, order);

		// It is possible that the buddy of source block is also available. Prepare to check and merge.
		PageDescriptor* to_merge = pgd;

		// Iteratively check up to MAX_ORDER-1, since on top order level you can't commit merge.
		for(int i = order; i < MAX_ORDER-1; i++)
		{
			// Declare and assign the "possible" buddy and iter slot for traversal.
			PageDescriptor* buddy = buddy_of(to_merge, i);
			PageDescriptor* slot = _free_areas[i];
			bool merge_flag = false;

			// Check if its buddy in "local" order level
			while(slot){
				// When find the buddy, merge them and update to_merge for continuous merging.
				if(slot == buddy){
					to_merge = *merge_block(&to_merge, i);

					// Make sure the product of the merge is correctly aligned in target level.
					assert(is_correct_alignment_for_order(to_merge, i+1));

					merge_flag = true;
					break;
				}

				slot = slot->next_free;
			}

			// Once find and merge buddy in "local" order level, could directly pass to next order level.
			if (!merge_flag) {
				break; 
			}
		}
	}

	/**
	 * Reserves a specific page, so that it cannot be allocated.
	 * @param pgd The page descriptor of the page to reserve.
	 * @return Returns TRUE if the reservation was successful, FALSE otherwise.
	 */
	bool reserve_page(PageDescriptor *pgd)
	{
		/** 
		 * Execution sequence: locating the block pgd in -> 
		 *                     iteratively split it down to zero-order ->
		 * 				       check and remove among zero-order blocks to reserve.
		 */ 
				
		// Initialize objects for locating the target pgd.
		PageDescriptor* locate_block = NULL;
		int locate_order = -1;
		
		// Iteratively searching to find the block and order pgd stay in _free_areas.
		for(int i = 0; i < MAX_ORDER; ++i)
		{	
			// Cancel the loop if you find where pgd is in.
			if (locate_block != NULL) {
				break;
			}

			// In each "local" order level, traversal search each block.
			PageDescriptor* slot = _free_areas[i];
			while(slot){
				// Try to locate by range checking
				if ((slot <= pgd) && (pgd < (slot + pages_per_block(i)))) {
					locate_block = slot;
					locate_order = i;
					break;
				}else{
					slot = slot->next_free;
				}
			}
		}

		// When you find the block pgd in, iteratively split it down to zero order to reserve it.
		if (locate_block != NULL) {			
			while(locate_order > 0){
				// Keep tracking on two product blocks.
				PageDescriptor* origin = split_block(&locate_block, locate_order);
				PageDescriptor* buddy = buddy_of(origin, locate_order-1);

				// Make sure the products are correctly aligned.
				assert(is_correct_alignment_for_order(origin, locate_order-1));
				assert(is_correct_alignment_for_order(buddy, locate_order-1));

				// pgd could be in one of them. Locate and update locate_block.
				if (pgd >= buddy) {
					locate_block = buddy;
				}else{
					locate_block = origin;
				}

				locate_order -= 1;
			}
		}

		// When you finish the spliting, check pgd among the single page blocks in _free_areas[0] and remove to reserve.
		PageDescriptor* zero_order_slot = _free_areas[0];
		while (zero_order_slot) {
			if (pgd == zero_order_slot) {
				remove_block(pgd, 0);
				return true;
			}else{
				zero_order_slot = zero_order_slot->next_free;
			}
		}

		// If the function executes to here, it means that the function can't find pgd in zero-order blocks. Which is fail.
		return false;
	}
	
	/**
	 * Initialises the allocation algorithm.
	 * @return Returns TRUE if the algorithm was successfully initialised, FALSE otherwise.
	 */
	bool init(PageDescriptor *page_descriptors, uint64_t nr_page_descriptors) override
	{
		mm_log.messagef(LogLevel::DEBUG, "Buddy Allocator Initialising pd=%p, nr=0x%lx", page_descriptors, nr_page_descriptors);
		
		// diff is the reminder of total pgd number moduls MAX_ORDER block size. 
		// For these pgds, they should be allocated to other order level instead of MAX_ORDER-1. 
		uint64_t diff =  nr_page_descriptors % pages_per_block(MAX_ORDER-1);

		// Initialize necessary counters.
		uint64_t load_to_MAX = nr_page_descriptors - diff;
		uint64_t block_size = pages_per_block(MAX_ORDER-1);
		uint64_t block_number = load_to_MAX / block_size;
		uint64_t allocated = 0;

		// Allocate to MAX_ORDER level if possible.
		if (block_number != 0) {
			_free_areas[MAX_ORDER-1] = page_descriptors; 
			PageDescriptor* to_load = _free_areas[MAX_ORDER-1];
			allocated += block_size;
			
			// Iteratively load pgds.
			for(uint64_t i = 0; i < block_number-1; i++)
			{	
				to_load->next_free = to_load + block_size;
				to_load = to_load->next_free;
				allocated += block_size;
			}

			// After loading, clear next_free pointer for last element in linked list.
			to_load->next_free = NULL;
		}

		// If no reminder, it means that all pgds are allocated.
		if (diff == 0) {
			return true;
		}else{
			// Otherwise, allocate the reminders. Iteratively allocate downward in each level.
			for(int i = MAX_ORDER-2; i >= 0; i--)
			{	
				// Once allocated all reminders, cancel the loop.
				if (diff == 0) {
					break;
				}
				
				// If current level perfectly fit the reminders, directly allocate and cancel.
				if (diff == pages_per_block(i)) {
					_free_areas[i] = page_descriptors + allocated;
					break;
				}else if(diff > pages_per_block(i)){
					// If reminders are more than current level block size, could allocate to this level.
					// Initialize necessary counters.
					uint64_t blk_size = pages_per_block(i);
					uint64_t buffer = diff;

					// Update the reminders number.
					diff = diff % blk_size;
					uint64_t blocks_to_fill = (buffer - diff) / blk_size;
					_free_areas[i] = page_descriptors + allocated;
					PageDescriptor* start = _free_areas[i];
					allocated += blk_size;
					
					// Iteratively allocate to currentlevel.
					for(uint64_t j = 0; j < blocks_to_fill-1; j++)
					{
						start->next_free = start + blk_size;
						start = start->next_free;
						allocated += blk_size;
					}

					// Clear the last next_free pointer
					start->next_free = NULL;
				}

				// Otherwise, the reminders number is less than the block size in current order level. Unable to allocate.
				// Directly pass to next order level and examine.
			}
		}
		return true;
	}

	/**
	 * Returns the friendly name of the allocation algorithm, for debugging and selection purposes.
	 */
	const char* name() const override { return "buddy"; }
	
	/**
	 * Dumps out the current state of the buddy system
	 */
	void dump_state() const override
	{
		// Print out a header, so we can find the output in the logs.
		mm_log.messagef(LogLevel::DEBUG, "BUDDY STATE:");
		
		// Iterate over each free area.
		for (unsigned int i = 0; i < ARRAY_SIZE(_free_areas); i++) {
			char buffer[256];
			snprintf(buffer, sizeof(buffer), "[%d] ", i);
						
			// Iterate over each block in the free area.
			PageDescriptor *pg = _free_areas[i];
			while (pg) {
				// Append the PFN of the free block to the output buffer.
				snprintf(buffer, sizeof(buffer), "%s%lx ", buffer, sys.mm().pgalloc().pgd_to_pfn(pg));
				pg = pg->next_free;
			}
			
			mm_log.messagef(LogLevel::DEBUG, "%s", buffer);
		}
	}

	
private:
	PageDescriptor *_free_areas[MAX_ORDER];
};

/* --- DO NOT CHANGE ANYTHING BELOW THIS LINE --- */

/*
 * Allocation algorithm registration framework
 */
RegisterPageAllocator(BuddyPageAllocator);
